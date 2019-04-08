/*
------------------------------------------------------------------

This file is part of a plugin for the Open Ephys GUI
Copyright (C) 2017 Translational NeuroEngineering Laboratory, MGH

------------------------------------------------------------------

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/

#include <cfloat>  // DBL_MAX
#include <cmath>   // sqrt
#include <cstring> // memcpy, memmove

#include "PhaseCalculator.h"
#include "PhaseCalculatorEditor.h"

namespace PhaseCalculator
{
    // ------- constants ------------

    // default passband width if pushing lowCut down or highCut up to fix invalid range,
    // and also the minimum for lowCut.
    // (meant to be larger than actual minimum floating-point eps)
    static const float passbandEps = 0.01F;

    // priority of the AR model calculating thread (0 = lowest, 10 = highest)
    static const int arPriority = 3;

    // "glitch limit" (how long of a segment is allowed to be unwrapped or smoothed, in samples)
    static const int glitchLimit = 200;

    // process length for real-time visualization ground truth hilbert transform
    // based on evaluation of phase error compared to offline processing
    static const int visHilbertLengthMs = 1024;
    static const int visMinDelayMs = 675;
    static const int visMaxDelayMs = 1000;


    /*** ReverseStack ***/
    ReverseStack::ReverseStack(int size)
        : freeSpace     (size)
        , headOffset    (size - 1)
    {
        resize(size);
    }

    void ReverseStack::reset()
    {
        const ScopedLock dataLock(getLock());
        freeSpace = size();
        headOffset = freeSpace - 1;
    }

    void ReverseStack::resetAndResize(int newSize)
    {
        const ScopedLock dataLock(getLock());
        resize(newSize);
        freeSpace = newSize;
        headOffset = newSize - 1;
    }

    bool ReverseStack::isFull() const
    {
        return freeSpace == 0;
    }

    void ReverseStack::enqueue(const float* source, int n)
    {
        const ScopedLock dataLock(getLock());

        // skip samples that can't be written
        int length = size();
        int nToSkip = jmax(0, n - length);
        int nToAdd = n - nToSkip;
        source += nToSkip;

        double* data = getRawDataPointer();

        for (int nLeft = nToAdd; nLeft > 0; --nLeft)
        {
            data[headOffset] = double(*(source++));
            headOffset = (headOffset ? headOffset : length) - 1;
        }

        freeSpace = jmax(0, freeSpace - nToAdd);
    }

    int ReverseStack::getHeadOffset() const
    {
        return headOffset;
    }

    void ReverseStack::unwrapAndCopy(double* dest, bool useLock) const
    {
        ScopedPointer<ScopedLock> lock;
        if (useLock)
        {
            lock = new ScopedLock(getLock());
        }

        int length = size();

        const double* block2Start = begin();
        int block2Size = (headOffset + 1) % length;

        const double* block1Start = block2Start + block2Size;
        int block1Size = length - block2Size;

        std::memcpy(dest, block1Start, block1Size * sizeof(double));
        std::memcpy(dest + block1Size, block2Start, block2Size * sizeof(double));
    }


    /**** channel info *****/
    ActiveChannelInfo::ActiveChannelInfo(const ChannelInfo& cInfo)
        : chanInfo(cInfo)
    {
        update();
    }

    void ActiveChannelInfo::update()
    {
        const Node& p = chanInfo.owner;
        int arOrder = p.getAROrder();
        float highCut = p.getHighCut();
        float lowCut = p.getLowCut();
        Band band = p.getBand();

        // update length of history based on sample rate
        // the history buffer should have enough samples to calculate phases for the viusalizer
        // with the proper Hilbert transform length AND train an AR model of the requested order,
        // using at least 1 second of data
        int newHistorySize = chanInfo.dsFactor * jmax(
            visHilbertLengthMs * Hilbert::fs / 1000,
            arOrder + 1,
            1 * Hilbert::fs);

        history.resetAndResize(newHistorySize);

        // set filter parameters
        Dsp::Params params;
        params[0] = chanInfo.sampleRate; // sample rate
        params[1] = 2;                          // order
        params[2] = (highCut + lowCut) / 2;     // center frequency
        params[3] = highCut - lowCut;           // bandwidth

        filter.setParams(params);

        arModeler.setParams(arOrder, newHistorySize, chanInfo.dsFactor);

        htState.resize(Hilbert::delay[band] * 2 + 1);

        // visualization stuff
        hilbertLengthMultiplier = Hilbert::fs * chanInfo.dsFactor / 1000;
        int visHilbertLength = visHilbertLengthMs * hilbertLengthMultiplier;

        if (visHilbertBuffer.getLength() != visHilbertLength)
        {
            visHilbertBuffer.resize(visHilbertLength);
            visForwardPlan = new FFTWPlan(visHilbertLength, &visHilbertBuffer, FFTW_MEASURE);
            visBackwardPlan = new FFTWPlan(visHilbertLength, &visHilbertBuffer, FFTW_BACKWARD, FFTW_MEASURE);
        }

        reverseFilter.setParams(params);

        reset();
    }

    void ActiveChannelInfo::reset()
    {
        history.reset();
        filter.reset();
        arModeler.reset();
        FloatVectorOperations::clear(htState.begin(), htState.size());
        dsOffset = chanInfo.dsFactor;
        lastComputedSample = 0;
        lastPhase = 0;
    }


    ChannelInfo::ChannelInfo(const Node& pc, int i)
        : chan          (i)
        , sampleRate    (0)
        , dsFactor      (0)
        , acInfo        (nullptr)
        , owner         (pc)
    {
        update();
    }

    void ChannelInfo::update()
    {
        const DataChannel* dataChannel = owner.getDataChannel(chan);
        if (dataChannel == nullptr)
        {
            jassertfalse;
            return;
        }

        sampleRate = dataChannel->getSampleRate();

        float fsMult = sampleRate / Hilbert::fs;
        float fsMultRound = std::round(fsMult);
        if (std::abs(fsMult - fsMultRound) < FLT_EPSILON)
        {
            // can be active - sample rate is multiple of Hilbert Fs
            dsFactor = int(fsMultRound);

            if (isActive())
            {
                acInfo->update();
            }
        }
        else
        {
            dsFactor = 0;
            deactivate(); // this channel can no longer be active.
        }
    }

    bool ChannelInfo::activate()
    {
        if (!isActive() && dsFactor != 0)
        {
            acInfo = new ActiveChannelInfo(*this);
        }

        return isActive();
    }

    void ChannelInfo::deactivate()
    {
        acInfo = nullptr;
    }

    bool ChannelInfo::isActive() const
    {
        return acInfo != nullptr;
    }

    /**** phase calculator node ****/
    Node::Node()
        : GenericProcessor("Phase Calculator")
        , Thread("AR Modeler")
        , calcInterval(50)
        , arOrder(20)
        , outputMode(PH)
        , visEventChannel(-1)
        , visContinuousChannel(-1)
    {
        setProcessorType(PROCESSOR_TYPE_FILTER);
        setBand(ALPHA_THETA, true);
    }

    Node::~Node() {}

    bool Node::hasEditor() const
    {
        return true;
    }


    AudioProcessorEditor* Node::createEditor()
    {
        editor = new Editor(this);
        return editor;
    }

    void Node::createEventChannels()
    {
        const DataChannel* visChannel = getDataChannel(visContinuousChannel);

        if (!visChannel)
        {
            visPhaseChannel = nullptr;
            return;
        }

        float sampleRate = visChannel->getSampleRate();

        EventChannel* chan = new EventChannel(EventChannel::DOUBLE_ARRAY, 1, 1, sampleRate, this);
        chan->setName(chan->getName() + ": PC visualized phase (deg.)");
        chan->setDescription("The accurate phase in degrees of each visualized event");
        chan->setIdentifier("phasecalc.visphase");

        // metadata storing source data channel
        MetaDataDescriptor sourceChanDesc(MetaDataDescriptor::UINT16, 3, "Source Channel",
            "Index at its source, Source processor ID and Sub Processor index of the channel that triggers this event",
            "source.channel.identifier.full");
        MetaDataValue sourceChanVal(sourceChanDesc);
        uint16 sourceInfo[3];
        sourceInfo[0] = visChannel->getSourceIndex();
        sourceInfo[1] = visChannel->getSourceNodeID();
        sourceInfo[2] = visChannel->getSubProcessorIdx();
        sourceChanVal.setValue(static_cast<const uint16*>(sourceInfo));
        chan->addMetaData(sourceChanDesc, sourceChanVal);

        visPhaseChannel = eventChannelArray.add(chan);
    }

    void Node::setParameter(int parameterIndex, float newValue)
    {
        switch (parameterIndex) {
        case RECALC_INTERVAL:
            calcInterval = int(newValue);
            break;

        case AR_ORDER:
            arOrder = int(newValue);
            updateActiveChannels();
            break;

        case BAND:
            setBand(Band(int(newValue)));
            break;

        case LOWCUT:
            setLowCut(newValue);
            break;

        case HIGHCUT:
            setHighCut(newValue);
            break;

        case OUTPUT_MODE:
        {
            OutputMode oldMode = outputMode;
            outputMode = OutputMode(int(newValue));
            if (oldMode == PH_AND_MAG || outputMode == PH_AND_MAG)
            {
                CoreServices::updateSignalChain(editor);  // add or remove channels if necessary
            }
            break;
        }

        case VIS_E_CHAN:
            jassert(newValue >= -1);
            visEventChannel = int(newValue);
            break;

        case VIS_C_CHAN:
            setVisContChan(int(newValue));
            break;
        }
    }

    void Node::process(AudioSampleBuffer& buffer)
    {
        // handle subprocessors, if any
        HashMap<int, uint16>::Iterator subProcIt(subProcessorMap);
        while (subProcIt.next())
        {
            uint32 fullSourceID = uint32(subProcIt.getKey());
            int subProcessor = subProcIt.getValue();
            uint64 sourceTimestamp = getSourceTimestamp(fullSourceID);
            uint32 sourceSamples = getNumSourceSamples(fullSourceID);
            setTimestampAndSamples(sourceTimestamp, sourceSamples, subProcessor);
        }

        // check for events to visualize
        bool hasCanvas = static_cast<Editor*>(getEditor())->canvas != nullptr;
        if (hasCanvas && visEventChannel > -1)
        {
            checkForEvents();
        }

        // iterate over active input channels
        Array<int> activeChans = getActiveInputs();
        int numActiveChans = activeChans.size();
        for (int ac = 0; ac < numActiveChans; ++ac)
        {
            ChannelInfo* chanInfo = channelInfo[activeChans[ac]];
            ActiveChannelInfo* acInfo = chanInfo->acInfo;

            int chan = chanInfo->chan;
            int nSamples = getNumSamples(chan);
            if (nSamples == 0) // nothing to do
            {
                continue;
            }

            // filter the data
            float* const wpIn = buffer.getWritePointer(chan);
            acInfo->filter.process(nSamples, &wpIn);

            // enqueue as much new data as can fit into history
            acInfo->history.enqueue(wpIn, nSamples);

            // calc phase and write out (only if AR model has been calculated)
            if (acInfo->history.isFull() && acInfo->arModeler.hasBeenFit())
            {
                // read current AR parameters safely (uses lock internally)
                acInfo->arModeler.getModel(localARParams);

                // use AR model to fill predSamps (which is downsampled) based on past data.
                int htDelay = Hilbert::delay[band];
                int stride = acInfo->chanInfo.dsFactor;

                double* pPredSamps = predSamps.getRawDataPointer();
                const double* pLocalParam = localARParams.getRawDataPointer();
                arPredict(acInfo->history, acInfo->dsOffset, pPredSamps, pLocalParam, htDelay + 1, stride, arOrder);

                // identify indices of current buffer to execute HT
                htInds.clearQuick();
                for (int i = stride - acInfo->dsOffset; i < nSamples; i += stride)
                {
                    htInds.add(i);
                }

                int htOutputSamps = htInds.size() + 1;
                if (htOutput.size() < htOutputSamps)
                {
                    htOutput.resize(htOutputSamps);
                }

                // execute tranformer on current buffer
                int kOut = -htDelay;
                for (int kIn = 0; kIn < htInds.size(); ++kIn, ++kOut)
                {
                    double samp = htFilterSamp(wpIn[htInds[kIn]], band, acInfo->htState);
                    if (kOut >= 0)
                    {
                        double rc = wpIn[htInds[kOut]];
                        double ic = htScaleFactor * samp;
                        htOutput.set(kOut, std::complex<double>(rc, ic));
                    }
                }

                // copy state to transform prediction without changing the end-of-buffer state
                htTempState = acInfo->htState;

                // execute transformer on prediction
                for (int i = 0; i <= htDelay; ++i, ++kOut)
                {
                    double samp = htFilterSamp(predSamps[i], band, htTempState);
                    if (kOut >= 0)
                    {
                        double rc = i == htDelay ? predSamps[0] : wpIn[htInds[kOut]];
                        double ic = htScaleFactor * samp;
                        htOutput.set(kOut, std::complex<double>(rc, ic));
                    }
                }

                // output with upsampling (interpolation)
                float* wpOut = buffer.getWritePointer(chan);
                float* wpOut2;
                if (outputMode == PH_AND_MAG)
                {
                    // second output channel
                    int outChan2 = getNumInputs() + ac;
                    jassert(outChan2 < buffer.getNumChannels());
                    wpOut2 = buffer.getWritePointer(outChan2);
                }

                kOut = 0;
                std::complex<double> prevCS = acInfo->lastComputedSample;
                std::complex<double> nextCS = htOutput[kOut];
                double prevPhase, nextPhase, phaseSpan, thisPhase;
                double prevMag, nextMag, magSpan, thisMag;
                bool needPhase = outputMode != MAG;
                bool needMag = outputMode != PH;

                if (needPhase)
                {
                    prevPhase = std::arg(prevCS);
                    nextPhase = std::arg(nextCS);
                    phaseSpan = circDist(nextPhase, prevPhase, Dsp::doublePi);
                }
                if (needMag)
                {
                    prevMag = std::abs(prevCS);
                    nextMag = std::abs(nextCS);
                    magSpan = nextMag - prevMag;
                }
                int subSample = acInfo->dsOffset % stride;

                for (int i = 0; i < nSamples; ++i, subSample = (subSample + 1) % stride)
                {
                    if (subSample == 0)
                    {
                        // update interpolation frame
                        prevCS = nextCS;
                        nextCS = htOutput[++kOut];

                        if (needPhase)
                        {
                            prevPhase = nextPhase;
                            nextPhase = std::arg(nextCS);
                            phaseSpan = circDist(nextPhase, prevPhase, Dsp::doublePi);
                        }
                        if (needMag)
                        {
                            prevMag = nextMag;
                            nextMag = std::abs(nextCS);
                            magSpan = nextMag - prevMag;
                        }
                    }

                    if (needPhase)
                    {
                        thisPhase = prevPhase + phaseSpan * subSample / stride;
                        thisPhase = circDist(thisPhase, 0, Dsp::doublePi);
                    }
                    if (needMag)
                    {
                        thisMag = prevMag + magSpan * subSample / stride;
                    }

                    switch (outputMode)
                    {
                    case MAG:
                        wpOut[i] = static_cast<float>(thisMag);
                        break;

                    case PH_AND_MAG:
                        wpOut2[i] = static_cast<float>(thisMag);
                        // fall through
                    case PH:
                        // output in degrees
                        wpOut[i] = static_cast<float>(thisPhase * (180.0 / Dsp::doublePi));
                        break;

                    case IM:
                        wpOut[i] = static_cast<float>(thisMag * std::sin(thisPhase));
                        break;
                    }
                }
                acInfo->lastComputedSample = prevCS;
                acInfo->dsOffset = ((acInfo->dsOffset + nSamples - 1) % stride) + 1;

                // unwrapping / smoothing
                if (outputMode == PH || outputMode == PH_AND_MAG)
                {
                    unwrapBuffer(wpOut, nSamples, acInfo->lastPhase);
                    smoothBuffer(wpOut, nSamples, acInfo->lastPhase);
                    acInfo->lastPhase = wpOut[nSamples - 1];
                }
            }
            else // fifo not full or AR model not ready
            {
                // just output zeros
                buffer.clear(chan, 0, nSamples);
            }

            // if this is the monitored channel for events, check whether we can add a new phase
            if (hasCanvas && chan == visContinuousChannel && acInfo->history.isFull())
            {
                calcVisPhases(acInfo, getTimestamp(chan) + getNumSamples(chan));
            }
        }
    }

    // starts thread when acquisition begins
    bool Node::enable()
    {
        if (isEnabled)
        {
            startThread(arPriority);

            // have to manually enable editor, I guess...
            Editor* editor = static_cast<Editor*>(getEditor());
            editor->enable();
        }

        return isEnabled;
    }

    bool Node::disable()
    {
        Editor* editor = static_cast<Editor*>(getEditor());
        editor->disable();

        signalThreadShouldExit();

        // reset states of active inputs
        for (auto chanInfo : channelInfo)
        {
            if (chanInfo->isActive())
            {
                chanInfo->acInfo->reset();
            }
        }

        // clear timestamp and phase queues
        while (!visTsBuffer.empty())
        {
            visTsBuffer.pop();
        }

        ScopedLock phaseLock(visPhaseBufferCS);
        while (!visPhaseBuffer.empty())
        {
            visPhaseBuffer.pop();
        }

        return true;
    }

    // thread routine
    void Node::run()
    {
        // collect enabled active channels and find maximum history length
        Array<ActiveChannelInfo*> activeChans;
        int maxHistoryLength = 0;
        for (auto chanInfo : channelInfo)
        {
            if (chanInfo->isActive())
            {
                activeChans.add(chanInfo->acInfo);
                maxHistoryLength = jmax(maxHistoryLength, chanInfo->acInfo->history.size());
            }
        }

        Array<double> reverseData;
        reverseData.resize(maxHistoryLength);

        uint32 startTime, endTime;
        while (!threadShouldExit())
        {
            startTime = Time::getMillisecondCounter();

            for (auto acInfo : activeChans)
            {
                if (!acInfo->history.isFull())
                {
                    continue;
                }

                // unwrap reversed history and add to temporary data array
                double* dataPtr = reverseData.getRawDataPointer();
                acInfo->history.unwrapAndCopy(dataPtr, true);

                // calculate parameters
                acInfo->arModeler.fitModel(reverseData);
            }

            endTime = Time::getMillisecondCounter();
            int remainingInterval = calcInterval - (endTime - startTime);
            if (remainingInterval >= 10) // avoid WaitForSingleObject
            {
                sleep(remainingInterval);
            }
        }
    }

    void Node::updateSettings()
    {
        // update arrays that store one entry per input
        int numInputs = getNumInputs();
        int prevNumInputs = channelInfo.size();

        int nToRemove = jmax(prevNumInputs - numInputs, 0);
        channelInfo.removeLast(nToRemove);

        updateAllChannels();

        for (int i = prevNumInputs; i < numInputs; ++i)
        {
            channelInfo.add(new ChannelInfo(*this, i));
        }

        // create new data channels if necessary
        updateSubProcessorMap();
        updateExtraChannels();

        if (outputMode == PH_AND_MAG)
        {
            // keep previously selected input channels from becoming selected extra channels
            deselectAllExtraChannels();
        }
    }


    Array<int> Node::getActiveInputs() const
    {
        Array<int> activeInputs;
        for (auto chanInfo : channelInfo)
        {
            if (chanInfo->isActive())
            {
                activeInputs.add(chanInfo->chan);
            }
        }
        return activeInputs;
    }


    bool Node::isGeneratesTimestamps() const
    {
        return true;
    }

    int Node::getNumSubProcessors() const
    {
        return subProcessorMap.size();
    }

    float Node::getSampleRate(int subProcessorIdx) const
    {
        jassert(subProcessorIdx < getNumSubProcessors());
        int chan = getDataChannelIndex(0, getNodeId(), subProcessorIdx);
        return getDataChannel(chan)->getSampleRate();
    }

    float Node::getBitVolts(int subProcessorIdx) const
    {
        jassert(subProcessorIdx < getNumSubProcessors());
        int chan = getDataChannelIndex(0, getNodeId(), subProcessorIdx);
        return getDataChannel(chan)->getBitVolts();
    }

    int Node::getFullSourceId(int chan)
    {
        const DataChannel* chanInfo = getDataChannel(chan);
        if (!chanInfo)
        {
            jassertfalse;
            return 0;
        }
        uint16 sourceNodeId = chanInfo->getSourceNodeID();
        uint16 subProcessorIdx = chanInfo->getSubProcessorIdx();
        return int(getProcessorFullId(sourceNodeId, subProcessorIdx));
    }

    int Node::getAROrder() const
    {
        return arOrder;
    }

    float Node::getHighCut() const
    {
        return highCut;
    }

    float Node::getLowCut() const
    {
        return lowCut;
    }

    Band Node::getBand() const
    {
        return band;
    }

    std::queue<double>& Node::getVisPhaseBuffer(ScopedPointer<ScopedLock>& lock)
    {
        lock = new ScopedLock(visPhaseBufferCS);
        return visPhaseBuffer;
    }

    void Node::saveCustomChannelParametersToXml(XmlElement* channelElement,
        int channelNumber, InfoObjectCommon::InfoObjectType channelType)
    {
        if (channelType == InfoObjectCommon::DATA_CHANNEL && channelNumber == visContinuousChannel)
        {
            channelElement->setAttribute("visualize", 1);
        }
    }

    void Node::loadCustomChannelParametersFromXml(XmlElement* channelElement,
        InfoObjectCommon::InfoObjectType channelType)
    {
        int chanNum = channelElement->getIntAttribute("number");

        if (chanNum < getNumInputs() && channelElement->hasAttribute("visualize"))
        {
            // The saved channel should be added to the dropdown at this point.
            setVisContChan(chanNum);
            static_cast<Editor*>(getEditor())->refreshVisContinuousChan();
        }
    }

    double Node::circDist(double x, double ref, double cutoff)
    {
        static const double twoPi = 2 * Dsp::doublePi;
        double xMod = std::fmod(x - ref, twoPi);
        double xPos = (xMod < 0 ? xMod + twoPi : xMod);
        return (xPos > cutoff ? xPos - twoPi : xPos);
    }

    // ------------ PRIVATE METHODS ---------------

    void Node::handleEvent(const EventChannel* eventInfo,
        const MidiMessage& event, int samplePosition)
    {
        if (visEventChannel < 0)
        {
            return;
        }

        if (Event::getEventType(event) == EventChannel::TTL)
        {
            TTLEventPtr ttl = TTLEvent::deserializeFromMessage(event, eventInfo);
            if (ttl->getChannel() == visEventChannel && ttl->getState())
            {
                // add timestamp to the queue for visualization
                juce::int64 ts = ttl->getTimestamp();
                jassert(visTsBuffer.empty() || visTsBuffer.back() <= ts);
                visTsBuffer.push(ts);
            }
        }
    }

    void Node::setBand(Band newBand, bool force)
    {
        if (!force && newBand == band) { return; }
        if (newBand < 0 || newBand >= NUM_BANDS)
        {
            jassertfalse;
            return;
        }

        band = newBand;

        // set low and high cut to the defaults for this band, making sure to notify the editor
        resetCutsToDefaults();

        // resize htState for each active channel, htTempState, and predSamps
        int delay = Hilbert::delay[band];
        htTempState.resize(delay * 2 + 1);
        predSamps.resize(delay + 1);

        updateActiveChannels();
    }

    void Node::resetCutsToDefaults()
    {
        const Array<float>& defaultBand = Hilbert::defaultBand[band];
        lowCut = defaultBand[0];
        highCut = defaultBand[1];

        auto editor = static_cast<Editor*>(getEditor());
        if (editor)
        {
            editor->refreshLowCut();
            editor->refreshHighCut();
        }

        updateScaleFactor();
        updateActiveChannels();
    }

    void Node::setLowCut(float newLowCut)
    {
        if (newLowCut == lowCut) { return; }

        auto editor = static_cast<Editor*>(getEditor());
        const Array<float>& validBand = Hilbert::validBand[band];

        if (newLowCut < validBand[0] || newLowCut >= validBand[1])
        {
            // invalid; don't set parameter and reset editor
            editor->refreshLowCut();
            CoreServices::sendStatusMessage("Low cut outside valid band of selected filter.");
            return;
        }

        lowCut = newLowCut;
        if (lowCut >= highCut)
        {
            // push highCut up
            highCut = jmin(lowCut + passbandEps, validBand[1]);
            editor->refreshHighCut();
        }

        updateScaleFactor();
        updateActiveChannels();
    }

    void Node::setHighCut(float newHighCut)
    {
        if (newHighCut == highCut) { return; }

        auto editor = static_cast<Editor*>(getEditor());
        const Array<float>& validBand = Hilbert::validBand[band];

        if (newHighCut <= validBand[0] || newHighCut > validBand[1])
        {
            // invalid; don't set parameter and reset editor
            editor->refreshHighCut();
            CoreServices::sendStatusMessage("High cut outside valid band of selected filter.");
            return;
        }

        highCut = newHighCut;
        if (highCut <= lowCut)
        {
            // push lowCut down
            lowCut = jmax(highCut - passbandEps, validBand[0]);
            editor->refreshLowCut();
        }

        updateScaleFactor();
        updateActiveChannels();
    }

    void Node::setVisContChan(int newChan)
    {
        if (newChan >= 0)
        {
            jassert(newChan < channelInfo.size() && channelInfo[newChan]->isActive());

            // disable event receival temporarily so we can flush the buffer
            int tempVisEventChan = visEventChannel;
            visEventChannel = -1;

            // clear timestamp queue
            while (!visTsBuffer.empty())
            {
                visTsBuffer.pop();
            }

            visEventChannel = tempVisEventChan;
        }
        
        visContinuousChannel = newChan;

        // If acquisition is stopped (and thus the new channel might be from a different subprocessor),
        // update signal chain. Sinks such as LFP Viewer should receive this information.
        if (!CoreServices::getAcquisitionStatus())
        {
            CoreServices::updateSignalChain(getEditor());
        }
    }

    void Node::updateScaleFactor()
    {
        htScaleFactor = getScaleFactor(band, lowCut, highCut);
    }

    void Node::unwrapBuffer(float* wp, int nSamples, float lastPhase)
    {
        for (int startInd = 0; startInd < nSamples - 1; startInd++)
        {
            float diff = wp[startInd] - (startInd == 0 ? lastPhase : wp[startInd - 1]);
            if (abs(diff) > 180)
            {
                // search forward for a jump in the opposite direction
                int endInd;
                int maxInd;
                if (diff < 0)
                    // for downward jumps, unwrap if there's a jump back up within glitchLimit samples
                {
                    endInd = -1;
                    maxInd = jmin(startInd + glitchLimit, nSamples - 1);
                }
                else
                    // for upward jumps, default to unwrapping until the end of the buffer, but stop if there's a jump back down sooner.
                {
                    endInd = nSamples;
                    maxInd = nSamples - 1;
                }
                for (int currInd = startInd + 1; currInd <= maxInd; currInd++)
                {
                    float diff2 = wp[currInd] - wp[currInd - 1];
                    if (abs(diff2) > 180 && ((diff > 0) != (diff2 > 0)))
                    {
                        endInd = currInd;
                        break;
                    }
                }

                // unwrap [startInd, endInd)
                for (int i = startInd; i < endInd; i++)
                {
                    wp[i] -= 360 * (diff / abs(diff));
                }

                if (endInd > -1)
                {
                    // skip to the end of this unwrapped section
                    startInd = endInd;
                }
            }
        }
    }

    void Node::smoothBuffer(float* wp, int nSamples, float lastPhase)
    {
        int actualGL = jmin(glitchLimit, nSamples - 1);
        float diff = wp[0] - lastPhase;
        if (diff < 0 && diff > -180)
        {
            // identify whether signal exceeds last sample of the previous buffer within glitchLimit samples.
            int endIndex = -1;
            for (int i = 1; i <= actualGL; i++)
            {
                if (wp[i] > lastPhase)
                {
                    endIndex = i;
                    break;
                }
                // corner case where signal wraps before it exceeds lastSample
                else if (wp[i] - wp[i - 1] < -180 && (wp[i] + 360) > lastPhase)
                {
                    wp[i] += 360;
                    endIndex = i;
                    break;
                }
            }

            if (endIndex != -1)
            {
                // interpolate points from buffer start to endIndex
                float slope = (wp[endIndex] - lastPhase) / (endIndex + 1);
                for (int i = 0; i < endIndex; i++)
                {
                    wp[i] = lastPhase + (i + 1) * slope;
                }
            }
        }
    }

    void Node::updateSubProcessorMap()
    {
        if (outputMode != PH_AND_MAG)
        {
            subProcessorMap.clear();
            return;
        }

        // fill map according to selected channels, and remove outdated entries.
        uint16 maxUsedIdx = 0;
        SortedSet<int> foundFullIds;
        Array<int> unmappedFullIds;

        Array<int> activeInputs = getActiveInputs();
        for (int chan : activeInputs)
        {
            const DataChannel* chanInfo = getDataChannel(chan);
            uint16 sourceNodeId = chanInfo->getSourceNodeID();
            uint16 subProcessorIdx = chanInfo->getSubProcessorIdx();
            int procFullId = int(getProcessorFullId(sourceNodeId, subProcessorIdx));
            foundFullIds.add(procFullId);

            if (subProcessorMap.contains(procFullId))
            {
                maxUsedIdx = jmax(maxUsedIdx, subProcessorMap[subProcessorIdx]);
            }
            else // add new entry for this source subprocessor
            {
                // try to match index if possible
                if (!subProcessorMap.containsValue(subProcessorIdx))
                {
                    subProcessorMap.set(procFullId, subProcessorIdx);
                    maxUsedIdx = jmax(maxUsedIdx, subProcessorIdx);
                }
                else
                {
                    unmappedFullIds.add(procFullId);
                }
            }
        }
        // assign remaining unmapped ids
        for (int id : unmappedFullIds)
        {
            subProcessorMap.set(id, ++maxUsedIdx);
        }

        // remove outdated entries
        Array<int> outdatedFullIds;
        HashMap<int, juce::uint16>::Iterator it(subProcessorMap);
        while (it.next())
        {
            int key = it.getKey();
            if (!foundFullIds.contains(key))
            {
                outdatedFullIds.add(key);
            }
        }
        for (int id : outdatedFullIds)
        {
            subProcessorMap.remove(id);
        }
    }

    void Node::updateExtraChannels()
    {
        // reset dataChannelArray to # of inputs
        int numInputs = getNumInputs();
        int numChannels = dataChannelArray.size();
        jassert(numChannels >= numInputs);
        dataChannelArray.removeLast(numChannels - numInputs);

        if (outputMode == PH_AND_MAG)
        {
            Array<int> activeInputs = getActiveInputs();
            for (int chan : activeInputs)
            {
                // see GenericProcessor::createDataChannelsByType
                DataChannel* baseChan = dataChannelArray[chan];
                int baseFullId = getFullSourceId(chan);

                DataChannel* newChan = new DataChannel(
                    baseChan->getChannelType(),
                    baseChan->getSampleRate(),
                    this,
                    subProcessorMap[baseFullId]);

                // rename to match base channel (implies that it contains magnitude data)
                newChan->setName(baseChan->getName() + "MAG");
                newChan->setBitVolts(baseChan->getBitVolts());
                newChan->addToHistoricString(getName());
                dataChannelArray.add(newChan);
            }
        }
        settings.numOutputs = dataChannelArray.size();
    }

    void Node::deselectChannel(int chan, bool warn)
    {
        jassert(chan >= 0 && chan < getTotalDataChannels());

        auto ed = getEditor();
        bool p, r, a;
        ed->getChannelSelectionState(chan, &p, &r, &a);
        ed->setChannelSelectionState(chan - 1, false, r, a);

        if (warn)
        {
            CoreServices::sendStatusMessage("Channel " + String(chan + 1) + " was deselected because" +
                " its sample rate is not a multiple of " + String(Hilbert::fs));
        }
    }

    void Node::deselectAllExtraChannels()
    {
        jassert(outputMode == PH_AND_MAG);
        Array<int> activeChans = getEditor()->getActiveChannels();
        int nInputs = getNumInputs();
        int nExtraChans = 0;
        for (int chan : activeChans)
        {
            if (chan < nInputs)
            {
                nExtraChans++;
            }
            else if (chan < nInputs + nExtraChans)
            {
                deselectChannel(chan, false);
            }
        }
    }

    void Node::calcVisPhases(ActiveChannelInfo* acInfo, juce::int64 sdbEndTs)
    {
        if (acInfo == nullptr)
        {
            jassertfalse;
            return;
        }

        int maxDelay = visMaxDelayMs * acInfo->hilbertLengthMultiplier;
        int minDelay = visMinDelayMs * acInfo->hilbertLengthMultiplier;
        int hilbertLength = visHilbertLengthMs * acInfo->hilbertLengthMultiplier;

        juce::int64 minTs = sdbEndTs - maxDelay;
        juce::int64 maxTs = sdbEndTs - minDelay;

        // discard any timestamps less than minTs
        while (!visTsBuffer.empty() && visTsBuffer.front() < minTs)
        {
            visTsBuffer.pop();
        }

        if (!visTsBuffer.empty() && visTsBuffer.front() <= maxTs)
        {
            // perform reverse filtering and Hilbert transform
            // don't need to use a lock here since it's the same thread as the one
            // that writes to it.
            double* wpHilbert = acInfo->visHilbertBuffer.getRealPointer();
            acInfo->history.unwrapAndCopy(wpHilbert, false);

            acInfo->reverseFilter.reset();
            acInfo->reverseFilter.process(hilbertLength, &wpHilbert);

            // un-reverse values
            acInfo->visHilbertBuffer.reverseReal(hilbertLength);

            acInfo->visForwardPlan->execute();
            hilbertManip(&acInfo->visHilbertBuffer, hilbertLength);
            acInfo->visBackwardPlan->execute();

            juce::int64 ts;
            ScopedLock phaseBufferLock(visPhaseBufferCS);
            while (!visTsBuffer.empty() && (ts = visTsBuffer.front()) <= maxTs)
            {
                visTsBuffer.pop();
                int delay = static_cast<int>(sdbEndTs - ts);
                std::complex<double> analyticPt = acInfo->visHilbertBuffer.getAsComplex(hilbertLength - delay);
                double phaseRad = std::arg(analyticPt);
                visPhaseBuffer.push(phaseRad);

                // add to event channel
                if (!visPhaseChannel)
                {
                    jassertfalse; // event channel should not be null here.
                    continue;
                }
                double eventData = phaseRad * 180.0 / Dsp::doublePi;
                juce::int64 eventTs = sdbEndTs - getNumSamples(acInfo->chanInfo.chan);
                BinaryEventPtr event = BinaryEvent::createBinaryEvent(visPhaseChannel, eventTs, &eventData, sizeof(double));
                addEvent(visPhaseChannel, event, 0);
            }
        }
    }

    void Node::updateAllChannels()
    {
        for (auto chanInfo : channelInfo)
        {
            bool wasActive = chanInfo->isActive();
            chanInfo->update();

            if (wasActive && !chanInfo->isActive())
            {
                // deselect if this channel just got deactivated
                deselectChannel(chanInfo->chan, true);
            }
        }
    }

    void Node::updateActiveChannels()
    {
        for (int ai : getActiveInputs())
        {
            jassert(channelInfo[ai] && channelInfo[ai]->isActive());
            channelInfo[ai]->acInfo->update();
        }
    }

    bool Node::activateInputChannel(int chan)
    {
        if (chan < 0 || chan >= channelInfo.size())
        {
            jassertfalse;
            return false;
        }

        jassert(!channelInfo[chan]->isActive()); // this shouldn't be called if it's already active.

        return channelInfo[chan]->activate();
    }

    void Node::deactivateInputChannel(int chan)
    {
        if (chan < 0 || chan >= channelInfo.size())
        {
            jassertfalse;
            return;
        }

        jassert(channelInfo.getUnchecked(chan)->isActive());
        channelInfo.getUnchecked(chan)->deactivate();
    }

    void Node::arPredict(const ReverseStack& history, int dsOffset, double* prediction,
        const double* params, int samps, int stride, int order)
    {
        const double* rpHistory = history.begin();
        int histSize = history.size();
        int histStart = history.getHeadOffset() + dsOffset;

        // s = index to write output
        for (int s = 0; s < samps; ++s)
        {
            prediction[s] = 0;

            // p = which AR param we are on
            for (int p = 0; p < order; ++p)
            {
                double pastSamp = p < s
                    ? prediction[s - 1 - p]
                    : rpHistory[(histStart + (p - s) * stride) % histSize];

                prediction[s] -= params[p] * pastSamp;
            }
        }
    }

    void Node::hilbertManip(FFTWArray* fftData, int n)
    {
        jassert(fftData->getLength() >= n);

        // Normalize DC and Nyquist, normalize and double positive freqs, and set negative freqs to 0.
        int lastPosFreq = (n + 1) / 2 - 1;
        int firstNegFreq = n / 2 + 1;
        int numPosNegFreqDoubles = lastPosFreq * 2; // sizeof(complex<double>) = 2 * sizeof(double)
        bool hasNyquist = (n % 2 == 0);

        std::complex<double>* wp = fftData->getComplexPointer();

        // normalize but don't double DC value
        wp[0] /= n;

        // normalize and double positive frequencies
        FloatVectorOperations::multiply(reinterpret_cast<double*>(wp + 1), 2.0 / n, numPosNegFreqDoubles);

        if (hasNyquist)
        {
            // normalize but don't double Nyquist frequency
            wp[lastPosFreq + 1] /= n;
        }

        // set negative frequencies to 0
        FloatVectorOperations::clear(reinterpret_cast<double*>(wp + firstNegFreq), numPosNegFreqDoubles);
    }

    double Node::getScaleFactor(Band band, double lowCut, double highCut)
    {
        double maxResponse = -DBL_MAX;
        double minResponse = DBL_MAX;

        Array<double> testFreqs({ lowCut, highCut });
        // also look at any magnitude response extrema that fall within the selected band
        for (double freq : Hilbert::extrema[band])
        {
            if (freq > lowCut && freq < highCut)
            {
                testFreqs.add(freq);
            }
        }

        // at each frequency, calculate the filter response
        int nCoefs = Hilbert::delay[band];
        for (double freq : testFreqs)
        {
            double normFreq = freq * Dsp::doublePi / (Hilbert::fs / 2);
            std::complex<double> response = 0;

            const double* transf = Hilbert::transformer[band].begin();
            for (int kCoef = 0; kCoef < nCoefs; ++kCoef)
            {
                double coef = transf[kCoef];

                // near component
                response += coef * std::polar(1.0, -(kCoef * normFreq));

                // mirrored component
                // there is no term for -nCoefs because that coefficient is 0.
                response -= coef * std::polar(1.0, -((2 * nCoefs - kCoef) * normFreq));
            }

            double absResponse = std::abs(response);
            maxResponse = jmax(maxResponse, absResponse);
            minResponse = jmin(minResponse, absResponse);
        }

        // scale factor is reciprocal of geometric mean of max and min
        return 1 / std::sqrt(minResponse * maxResponse);
    }

    double Node::htFilterSamp(double input, Band band, Array<double>& state)
    {
        double* state_p = state.getRawDataPointer();

        // initialize new state entry
        int nCoefs = Hilbert::delay[band];
        int order = nCoefs * 2;
        jassert(order == state.size() - 1);
        state_p[order] = 0;

        // incorporate new input
        const double* transf = Hilbert::transformer[band].begin();
        for (int kCoef = 0; kCoef < nCoefs; ++kCoef)
        {
            double val = input * transf[kCoef];
            state_p[kCoef] += val;          // near component
            state_p[order - kCoef] -= val;  // mirrored component
        }

        // output and shift state
        double sampOut = state_p[0];
        std::memmove(state_p, state_p + 1, order * sizeof(double));
        return sampOut;
    }
}


