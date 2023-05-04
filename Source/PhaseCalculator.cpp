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
            data[headOffset--] = double(*(source++));
            if (headOffset < 0)
            {
                headOffset = length - 1;
            }
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
    ActiveChannelInfo::ActiveChannelInfo(const ChannelInfo* cInfo)
        : chanInfo(cInfo)
    {

        bufferResizeThread = std::make_unique<BufferResizeThread>(&visHilbertBuffer);

        update();
    }

    void ActiveChannelInfo::update()
    {
        const DataStream* ds = chanInfo->stream;
        int arOrder = ds->getParameter("ar_order")->getValue();
        float highCut = ds->getParameter("high_cut")->getValue();
        float lowCut = ds->getParameter("low_cut")->getValue();
        Band band = (Band)static_cast<CategoricalParameter*>(ds->getParameter("freq_range"))->getSelectedIndex();

        // update length of history based on sample rate
        // the history buffer should have enough samples to calculate phases for the viusalizer
        // with the proper Hilbert transform length AND train an AR model of the requested order,
        // using at least 1 second of data
        int newHistorySize = chanInfo->dsFactor * jmax(
            visHilbertLengthMs * Hilbert::fs / 1000,
            arOrder + 1,
            1 * Hilbert::fs);

        LOGC("PhaseCalculator: Resetting history size");
        history.resetAndResize(newHistorySize);

        // set filter parameters
        for (auto filt : { &filter, &reverseFilter })
        {
            filt->setup(
                2,                      // order
                chanInfo->sampleRate,    // sample rate
                (highCut + lowCut) / 2, // center frequency
                highCut - lowCut);      // bandwidth
        }

        LOGC("PhaseCalculator: Setting filter parameters");
        arModeler.setParams(arOrder, newHistorySize, chanInfo->dsFactor);

        LOGC("PhaseCalculator: Resizing hilbert state");
        htState.resize(Hilbert::delay[band] * 2 + 1);

        // visualization stuff
        hilbertLengthMultiplier = Hilbert::fs * chanInfo->dsFactor / 1000;

        LOGC("PhaseCalculator: Resizing visualization buffer to ", visHilbertLengthMs * hilbertLengthMultiplier);

        if (bufferResizeThread->isThreadRunning())
        {
            bufferResizeThread->waitForThreadToExit(5000);
        }
        bufferResizeThread->setSize(visHilbertLengthMs * hilbertLengthMultiplier);
        bufferResizeThread->startThread();
        bufferResizeThread->waitForThreadToExit(10000);

        LOGC("PhaseCalculator: Resetting info");
        reset();
    }

    void ActiveChannelInfo::reset()
    {
        history.reset();
        filter.reset();
        arModeler.reset();
        FloatVectorOperations::clear(htState.begin(), htState.size());
        interpCountdown = 0;
        lastComputedPhase = 0;
        lastComputedMag = 0;
        lastPhase = 0;
    }


    ChannelInfo::ChannelInfo(const DataStream* ds, int i)
        : chan          (i)
        , sampleRate    (0)
        , dsFactor      (0)
        , isActivated   (false)
        , stream         (ds)
    {
        update();
    }

    void ChannelInfo::update()
    {
        ContinuousChannel* contChannel = stream->getContinuousChannels().getUnchecked(chan); 
        
        if (contChannel == nullptr)
        {
            jassertfalse;
            return;
        }

        sampleRate = contChannel->getSampleRate();

        float fsMult = sampleRate / Hilbert::fs;
        float fsMultRound = std::round(fsMult);
        
        if (std::abs(fsMult - fsMultRound) < FLT_EPSILON)
        {
            // can be active - sample rate is multiple of Hilbert Fs
            dsFactor = int(fsMultRound);

            if (isActivated)
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
        if (!isActivated && dsFactor != 0)
        {
            acInfo.reset(new ActiveChannelInfo(this));
            isActivated = true;
        }

        return isActivated;
    }

    void ChannelInfo::deactivate()
    {
        isActivated = false;
    }

    bool ChannelInfo::isActive() const
    {
        return isActivated;
    }


    /******** Phase Calculator Stream Settings *****/
    Settings::Settings():
        calcInterval(50),
        arOrder(20),
        visContinuousChannel(-1),
        visEventChannel(-1)
    {
        channelInfo.clear();
        setBand(Band(0), true);
    }

    Array<int> Settings::getActiveInputs() const
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

    void Settings::updateActiveChannels()
    {
        for (int ai : getActiveInputs())
        {
            jassert(channelInfo[ai] && channelInfo[ai]->isActive());
            channelInfo[ai]->acInfo->update();
        }
    }

    bool Settings::activateInputChannel(int chan)
    {
        if (chan < 0 || chan >= channelInfo.size())
        {
            jassertfalse;
            return false;
        }

        //jassert(!channelInfo[chan]->isActive()); // this shouldn't be called if it's already active.

        return channelInfo[chan]->activate();
    }

    void Settings::deactivateInputChannel(int chan)
    {
        if (chan < 0 || chan >= channelInfo.size())
        {
            jassertfalse;
            return;
        }

        //jassert(channelInfo.getUnchecked(chan)->isActive());
        channelInfo.getUnchecked(chan)->deactivate();
    }

    void Settings::setBand(Band newBand, bool force)
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

    }

    void Settings::resetCutsToDefaults()
    {
        const Array<float>& defaultBand = Hilbert::defaultBand[band];
        lowCut = defaultBand[0];
        highCut = defaultBand[1];

        updateScaleFactor();
        updateActiveChannels();
    }

    void Settings::updateScaleFactor()
    {
        htScaleFactor = getScaleFactor(band, lowCut, highCut);
    }

    double Settings::getScaleFactor(Band band, double lowCut, double highCut)
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


    /**** phase calculator node ****/
    Node::Node()
        : GenericProcessor("Phase Calculator")
        , Thread("AR Modeler")
    {
        setProcessorType(Plugin::Processor::FILTER);

        selectedStream = 0;
        activeChansNeedsUpdate = true;

        StringArray bands;
        for (int b = 0; b < NUM_BANDS - 1; ++b)
            bands.add(Hilbert::bandName[b]);

        addCategoricalParameter(Parameter::STREAM_SCOPE, 
                                "freq_range", 
                                "Each option corresponds internally to a Hilbert transformer that is optimized for this frequency range."
                                + String("After selecting a range, you can adjust ") +
                                "'low' and 'high' to filter to any passband within this range.",
                                bands, 0 );

        String desc;

        const Array<float>& defaultBand = Hilbert::defaultBand[0];
        
        addFloatParameter(Parameter::STREAM_SCOPE, "low_cut", "Low Cut", defaultBand[0], 0.0f, 1000.0f, 1.0f);
        
        addFloatParameter(Parameter::STREAM_SCOPE, "high_cut", "High Cut", defaultBand[1], 0.0f, 1000.0f, 1.0f);
        
        desc = "Time to wait between calls to update the autoregressive models"; 
        addIntParameter(Parameter::STREAM_SCOPE, "ar_refresh", desc, 50, 0, 10000);
        
        desc = "Order of the autoregressive models used to predict future data";
        addIntParameter(Parameter::STREAM_SCOPE, "ar_order", desc, 20, 1, 1000);
        
        addSelectedChannelsParameter(Parameter::STREAM_SCOPE, "Channels", "Selectable Channels", std::numeric_limits<int>::max());

        addIntParameter(Parameter::STREAM_SCOPE, "vis_cont", "Phase calculation channel", -1, -1, 1000);
        addIntParameter(Parameter::STREAM_SCOPE, "vis_event", "Event channel to plot phases", -1, -1, 1000);
    }


    AudioProcessorEditor* Node::createEditor()
    {
        editor = std::make_unique<Editor>(this);
        return editor.get();
    }


    void Node::process(AudioBuffer<float>& buffer)
    {

        // check for events to visualize
        bool hasCanvas = static_cast<Editor*>(getEditor())->canvas != nullptr;

        if (hasCanvas && settings[selectedStream]->visEventChannel > -1)
        {
            checkForEvents();
        }

        for (auto stream : dataStreams)
        {
            
            if ((*stream)["enable_stream"]
                && stream->getStreamId() == selectedStream)
            {
                // iterate over active input channels
                Array<int> activeChans = settings[stream->getStreamId()]->getActiveInputs();
                int numActiveChans = activeChans.size();

                int nSamples = getNumSamplesInBlock(stream->getStreamId());

                for (int ac = 0; ac < numActiveChans; ++ac)
                {
                    ChannelInfo* chanInfo = settings[stream->getStreamId()]->channelInfo[activeChans[ac]];
                    ActiveChannelInfo* acInfo = chanInfo->acInfo.get();

                    int chan = stream->getContinuousChannels().getUnchecked(chanInfo->chan)->getGlobalIndex();
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
                        int htDelay = Hilbert::delay[settings[stream->getStreamId()]->band];
                        int stride = chanInfo->dsFactor;

                        double* pPredSamps = settings[stream->getStreamId()]->predSamps.getRawDataPointer();
                        const double* pLocalParam = localARParams.getRawDataPointer();
                        arPredict(acInfo->history, acInfo->interpCountdown, pPredSamps, pLocalParam,
                            htDelay + 1, stride, settings[stream->getStreamId()]->arOrder);

                        // identify indices of current buffer to execute HT
                        htInds.clearQuick();
                        for (int i = acInfo->interpCountdown; i < nSamples; i += stride)
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
                            double samp = htFilterSamp(wpIn[htInds[kIn]], settings[stream->getStreamId()]->band, acInfo->htState);
                            if (kOut >= 0)
                            {
                                double rc = wpIn[htInds[kOut]];
                                double ic = settings[stream->getStreamId()]->htScaleFactor * samp;
                                htOutput.set(kOut, std::complex<double>(rc, ic));
                            }
                        }

                        // copy state to transform prediction without changing the end-of-buffer state
                        settings[stream->getStreamId()]->htTempState = acInfo->htState;

                        // execute transformer on prediction
                        for (int i = 0; i <= htDelay; ++i, ++kOut)
                        {
                            double samp = htFilterSamp(settings[stream->getStreamId()]->predSamps[i], 
                                            settings[stream->getStreamId()]->band,
                                            settings[stream->getStreamId()]->htTempState);
                            if (kOut >= 0)
                            {
                                double rc = i == htDelay ? settings[stream->getStreamId()]->predSamps[0] : wpIn[htInds[kOut]];
                                double ic = settings[stream->getStreamId()]->htScaleFactor * samp;
                                htOutput.set(kOut, std::complex<double>(rc, ic));
                            }
                        }

                        // output with upsampling (interpolation)
                        float* wpOut = buffer.getWritePointer(chan);

                        double nextComputedPhase, phaseStep;

                        nextComputedPhase = std::arg(htOutput[0]);
                        phaseStep = circDist(nextComputedPhase, acInfo->lastComputedPhase, Dsp::doublePi) / stride;

                        for (int i = 0, frame = 0; i < nSamples; ++i, --acInfo->interpCountdown)
                        {
                            if (acInfo->interpCountdown == 0)
                            {
                                // update interpolation frame
                                ++frame;
                                acInfo->interpCountdown = stride;

                                acInfo->lastComputedPhase = nextComputedPhase;
                                nextComputedPhase = std::arg(htOutput[frame]);

                            }

                            double thisPhase, thisMag;
                            thisPhase = circDist(nextComputedPhase, phaseStep * acInfo->interpCountdown, Dsp::doublePi);
                            wpOut[i] = float(thisPhase * (180.0 / Dsp::doublePi));

                        }

                        // unwrapping / smoothing
                        unwrapBuffer(wpOut, nSamples, acInfo->lastPhase);
                        smoothBuffer(wpOut, nSamples, acInfo->lastPhase);
                        acInfo->lastPhase = wpOut[nSamples - 1];
                    }
                    else // fifo not full or AR model not ready
                    {
                        // just output zeros
                        buffer.clear(chan, 0, nSamples);
                    }

                    // if this is the monitored channel for events, check whether we can add a new phase
                    if (hasCanvas
                        && chanInfo->chan == settings[stream->getStreamId()]->visContinuousChannel
                        && acInfo->history.isFull())
                    {
                        calcVisPhases(acInfo, getFirstSampleNumberForBlock(stream->getStreamId()));
                    }
                }
            }
        }
    }

    bool Node::startAcquisition()
    {
        if (isEnabled)
        {

            // wait for initialization threads to finish
            for (auto stream : getDataStreams())
            {
                for (auto chanInfo : settings[stream->getStreamId()]->channelInfo)
                {
                    if (chanInfo->isActive())
                    {
                        chanInfo->acInfo->waitForThreadToExit();
                    }
                }
            }

            activeChansNeedsUpdate = true;
            startThread(arPriority);

            // have to manually enable editor, I guess...
            Editor* editor = static_cast<Editor*>(getEditor());
            editor->enable();
        }

        return isEnabled;
    }

    bool Node::stopAcquisition()
    {
        Editor* editor = static_cast<Editor*>(getEditor());
        editor->disable();

        signalThreadShouldExit();

        // reset states of active inputs
        for(auto stream : getDataStreams())
        {
            for (auto chanInfo : settings[stream->getStreamId()]->channelInfo)
            {
                if (chanInfo->isActive())
                {
                    chanInfo->acInfo->reset();
                }
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

    void Node::setSelectedStream(uint16 streamID)
    {
        selectedStream = streamID;
        activeChansNeedsUpdate = true;
    }

    // thread routine
    void Node::run()
    {
        Array<ActiveChannelInfo*> activeChans;
        int maxHistoryLength = 0;

        Array<double> reverseData;

        uint32 startTime, endTime;
        while (!threadShouldExit())
        {
            startTime = Time::getMillisecondCounter();

            // collect enabled active channels and find maximum history length
            if(activeChansNeedsUpdate)
            {
                activeChans.clear();
                maxHistoryLength = 0;

                for (auto chanInfo : settings[selectedStream]->channelInfo)
                {
                    if (chanInfo->isActive())
                    {
                        activeChans.add(chanInfo->acInfo.get());
                        maxHistoryLength = jmax(maxHistoryLength, chanInfo->acInfo->history.size());
                    }
                }

                reverseData.resize(maxHistoryLength);
                activeChansNeedsUpdate = false;
            }

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
            int remainingInterval = settings[selectedStream]->calcInterval - (endTime - startTime);
            if (remainingInterval >= 10) // avoid WaitForSingleObject
            {
                sleep(remainingInterval);
            }
        }
    }

    void Node::updateSettings()
    {
        settings.update(getDataStreams());

        for (auto stream : getDataStreams())
        {
            settings[stream->getStreamId()]->channelInfo.clear();

            for (int i = 0; i < stream->getChannelCount(); ++i)
            {
                ChannelInfo* addChanInfo = new ChannelInfo(stream, i);
                settings[stream->getStreamId()]->channelInfo.add(addChanInfo);
            }

            parameterValueChanged(stream->getParameter("Channels"));
        }

    }


    void Node::parameterValueChanged(Parameter* param)
    {
        juce::uint16 paramStreamId = param->getStreamId();
        auto stream = getDataStream(paramStreamId);

        LOGD("[PhaseCalc] Parameter value changed ", paramStreamId, " : ", param->getName(), " : ", (int)param->getValue());

        if (param->getName().equalsIgnoreCase("Channels"))
        {   
            auto paramValue = static_cast<SelectedChannelsParameter*>(param)->getValue();

            LOGD("[PhaseCalc] Selected Channels Size: ", paramValue.getArray()->size());

            bool paramNeedsUpdate = false;

            for (int i = 0; i < getDataStream(paramStreamId)->getChannelCount(); i++)
            {
                // mark channel as activated or deactivated
                if (paramValue.getArray()->contains(i))
                {
                    // check whether channel can be activated
                    if (!settings[paramStreamId]->activateInputChannel(i))
                    {
                        LOGD("Failed to activate input channel ", i);
                        paramValue.getArray()->removeFirstMatchingValue(i);
                        paramNeedsUpdate = true;
                        continue;
                    }
                }
                else
                {
                    settings[paramStreamId]->deactivateInputChannel(i);
                }
            }

            setSelectedStream(paramStreamId);

            if(paramNeedsUpdate)
            {
                param->currentValue = paramValue;
                getEditor()->updateView();
            }

            getEditor()->updateVisualizer();
            
        }
        else if(param->getName().equalsIgnoreCase("freq_range"))
        {
            settings[paramStreamId]->setBand(Band((int)param->getValue()), true);

            if(stream->getParameter("low_cut") != nullptr && stream->getParameter("high_cut") != nullptr)
            {
                stream->getParameter("low_cut")->setNextValue(settings[paramStreamId]->lowCut);
                stream->getParameter("high_cut")->setNextValue(settings[paramStreamId]->highCut);
                settings[paramStreamId]->updateActiveChannels();
            }
        }
        else if(param->getName().equalsIgnoreCase("ar_refresh"))
        {
            settings[paramStreamId]->calcInterval = param->getValue();
        }
        else if(param->getName().equalsIgnoreCase("ar_order"))
        {
            settings[paramStreamId]->arOrder = param->getValue();
            settings[paramStreamId]->updateActiveChannels();
        }
        else if(param->getName().equalsIgnoreCase("low_cut"))
        {
            float newLowCut = (float)param->getValue();
            if (newLowCut == settings[paramStreamId]->lowCut) 
                return; 

            const Array<float>& validBand = Hilbert::validBand[settings[paramStreamId]->band];

            if (newLowCut < validBand[0] || newLowCut >= validBand[1])
            {
                // invalid; don't set parameter and reset editor
                CoreServices::sendStatusMessage("Low cut outside valid band of selected filter.");
                param->restorePreviousValue();
                return;
            }

            settings[paramStreamId]->lowCut = newLowCut;

            if (newLowCut >= settings[paramStreamId]->highCut)
            {
                // push highCut up
                stream->getParameter("high_cut")->setNextValue(jmin(newLowCut + passbandEps, validBand[1]));
            }
            else
            {
                settings[paramStreamId]->updateScaleFactor();
                settings[paramStreamId]->updateActiveChannels();
            }
        }
        else if(param->getName().equalsIgnoreCase("high_cut"))
        {
            float newHighCut = (float)param->getValue();
            if (newHighCut == settings[paramStreamId]->highCut)
                return;

            const Array<float>& validBand = Hilbert::validBand[settings[paramStreamId]->band];

            if (newHighCut <= validBand[0] || newHighCut > validBand[1])
            {

                // invalid; don't set parameter and reset editor
                CoreServices::sendStatusMessage("High cut outside valid band of selected filter.");
                param->restorePreviousValue();
                return;
            }

            settings[paramStreamId]->highCut = newHighCut;
            if (newHighCut <= settings[paramStreamId]->lowCut)
            {
                // push lowCut down
                stream->getParameter("low_cut")->setNextValue(jmax(newHighCut - passbandEps, validBand[0]));
            }
            else
            {
                settings[paramStreamId]->updateScaleFactor();
                settings[paramStreamId]-> updateActiveChannels();
            }
        }
        else if(param->getName().equalsIgnoreCase("vis_cont"))
        {
            setVisContChan((int)param->getValue());
        }
        else if(param->getName().equalsIgnoreCase("vis_event"))
        {
            jassert((int)param->getValue() >= -1);
            settings[paramStreamId]->visEventChannel = (int)param->getValue();
        }
        else
        {
            //do nothing
        }
    }

    bool Node::tryToReadVisPhases(std::queue<double>& other)
    {
        const ScopedTryLock lock(visPhaseBufferCS);
        if (!lock.isLocked())
        {
            return false;
        }

        while (!visPhaseBuffer.empty())
        {
            other.push(visPhaseBuffer.front());
            visPhaseBuffer.pop();
        }
        return true;
    }

    double Node::circDist(double x, double ref, double cutoff)
    {
        static const double twoPi = 2 * Dsp::doublePi;
        double xMod = std::fmod(x - ref, twoPi);
        double xPos = (xMod < 0 ? xMod + twoPi : xMod);
        return (xPos > cutoff ? xPos - twoPi : xPos);
    }

    // ------------ PRIVATE METHODS ---------------

    void Node::handleTTLEvent(TTLEventPtr event)
    {
        if (settings[selectedStream]->visEventChannel < 0)
        {
            return;
        }

        if (event->getEventType() == EventChannel::TTL)
        {
            if (event->getStreamId() == selectedStream
                && event->getLine() == settings[selectedStream]->visEventChannel
                && event->getState())
            {
                // add timestamp to the queue for visualization
                juce::int64 ts = event->getSampleNumber();
                jassert(visTsBuffer.empty() || visTsBuffer.back() <= ts);
                visTsBuffer.push(ts);
            }
        }
    }


    void Node::setVisContChan(int newChan)
    {
        if (newChan >= 0)
        {
            // jassert(newChan < channelInfo.size() && channelInfo[newChan]->isActive());

            // disable event receival temporarily so we can flush the buffer
            int tempVisEventChan = settings[selectedStream]->visEventChannel;
            settings[selectedStream]->visEventChannel = -1;

            // clear timestamp queue
            while (!visTsBuffer.empty())
            {
                visTsBuffer.pop();
            }

            settings[selectedStream]->visEventChannel = tempVisEventChan;
        }
        
        settings[selectedStream]->visContinuousChannel = newChan;
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


    Array<int> Node::getActiveChannels()
    {
        Array<int> activeInputs;

        if (selectedStream != 0)
        {
            for(int i : settings[selectedStream]->getActiveInputs())
                activeInputs.add(i);
        }
        
        return activeInputs;
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

            // Hilbert transform!
            acInfo->visHilbertBuffer.hilbert();

            juce::int64 ts;
            ScopedLock phaseBufferLock(visPhaseBufferCS);
            while (!visTsBuffer.empty() && (ts = visTsBuffer.front()) <= maxTs)
            {
                visTsBuffer.pop();
                int delay = static_cast<int>(sdbEndTs - ts);
                std::complex<double> analyticPt = acInfo->visHilbertBuffer.getAsComplex(hilbertLength - delay);
                double phaseRad = std::arg(analyticPt);
                visPhaseBuffer.push(phaseRad);

            }
        }
    }

    void Node::arPredict(const ReverseStack& history, int interpCountdown, double* prediction,
        const double* params, int samps, int stride, int order)
    {
        const double* rpHistory = history.begin();
        int histSize = history.size();
        int histStart = history.getHeadOffset() + stride - interpCountdown;

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


