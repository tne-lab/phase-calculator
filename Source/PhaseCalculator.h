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

#ifndef PHASE_CALCULATOR_H_INCLUDED
#define PHASE_CALCULATOR_H_INCLUDED

/*

The Phase Calculator filters selected input channels to the given passband,
generates an estimate of the phase of each band-limited signal in degrees
using the Hilbert transform, and outputs this as a continuous stream.
There are also options to output the magnitude or imaginary component of the analytic signal.

The visualizer panel serves a complementary role. When given a continuous channel
(that is enabled for processing, i.e. selected in the param window) and an event
channel to monitor, it calculates a delayed but precise estimate of the phase of the
filtered continuous channel at the time of each event onset and adds it to the rose plot.
If the events represent stimulation times, for example, this allows you to monitor the
accuracy of phase-locked stimulation in real time.

@see GenericProcessor, PhaseCalculatorEditor

*/

#include <ProcessorHeaders.h>
#include <DspLib.h>  // Filtering
#include <FFTWWrapper.h>   // Fourier transform

#include <queue>
#include <utility>  // pair

#include "ARModeler.h"     // Autoregressive modeling
#include "HTransformers.h" // Hilbert transformers & frequency bands

namespace PhaseCalculator
{
    // forward declarations
    struct ChannelInfo;
    class Node;

    // parameter indices
    enum Param
    {
        RECALC_INTERVAL,
        AR_ORDER,
        BAND,
        LOWCUT,
        HIGHCUT,
        OUTPUT_MODE,
        VIS_E_CHAN,
        VIS_C_CHAN
    };

    class ReverseStack : public Array<double, CriticalSection>
    {
    public:
        ReverseStack(int size = 0);

        // Just resets the free space
        void reset();

        // Resets free space and resizes the array. There's no way to resize
        // the array while keeping the data contained in it.
        void resetAndResize(int newSize);

        bool isFull() const;

        void enqueue(const float* source, int n);

        int getHeadOffset() const;

        // copies to given array, with first element corresponding to most recent data point.
        void unwrapAndCopy(double* dest, bool useLock) const;

    private:
        int freeSpace;
        int headOffset; // offset of write head from start of buffer

        JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(ReverseStack);
    };


    struct ActiveChannelInfo
    {
        ActiveChannelInfo(const ChannelInfo& cInfo);

        void update();

        // reset to perform after end of acquisition or update
        void reset();

        // filter design copied from FilterNode
        using BandpassFilter = Dsp::SimpleFilter
            <Dsp::Butterworth::BandPass // filter type
            <2>,                        // order
            1,                          // number of channels
            Dsp::DirectFormII>;         // realization

        ReverseStack history;

        BandpassFilter filter;

        ARModeler arModeler;

        Array<double> htState;

        // number of samples until a new non-interpolated output. e.g. if this
        // equals 1 after a buffer is processed, then there is one interpolated
        // sample in the next buffer, and then the second sample will be computed.
        // in range [0, dsFactor).
        int interpCountdown;

        // last non-interpolated ("computed") transformer output
        double lastComputedPhase;
        double lastComputedMag;

        // last phase output, for glitch correction
        float lastPhase;

        // for visualization:
        int hilbertLengthMultiplier;
        FFTWTransformableArray visHilbertBuffer;
        BandpassFilter reverseFilter;

        const ChannelInfo& chanInfo;

    private:
        JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(ActiveChannelInfo);
    };

    struct ChannelInfo
    {
        ChannelInfo(const Node& pc, int i);

        void update();

        // returns false on failure. does nothing if already active.
        bool activate();

        void deactivate();

        bool isActive() const;

        // index of channel among owner's data channels
        int chan;

        float sampleRate;

        // 0 if sample rate is not a multiple of Hilbert::FS (in this case it cannot be activated.)
        int dsFactor;

        // info for ongoing phase calculation - null if non-active.
        ScopedPointer<ActiveChannelInfo> acInfo;
        const Node& owner;

    private:
        JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(ChannelInfo);
    };

    // Output mode - corresponds to itemIDs on the ComboBox
    enum OutputMode { PH = 1, MAG, PH_AND_MAG, IM };

    class Node : public GenericProcessor, public Thread
    {
        friend class Editor;
    public:

        Node();
        ~Node();

        bool hasEditor() const override;

        AudioProcessorEditor* createEditor() override;

        void createEventChannels() override;

        void setParameter(int parameterIndex, float newValue) override;

        void process(AudioSampleBuffer& buffer) override;

        bool enable() override;
        bool disable() override;

        // thread code - recalculates AR parameters.
        void run() override;

        // handle changing number of channels
        void updateSettings() override;

        // Returns array of active channels that only includes inputs (not extra outputs)
        Array<int> getActiveInputs() const;

        // ----- to create new channels for multiple outputs -------
        bool isGeneratesTimestamps() const override;
        int getNumSubProcessors() const override;
        float getSampleRate(int subProcessorIdx = 0) const override;
        float getBitVolts(int subProcessorIdx = 0) const override;

        // miscellaneous helper
        int getFullSourceId(int chan);

        // getters
        int getAROrder() const;
        float getHighCut() const;
        float getLowCut() const;
        Band getBand() const;

        // reads from the visPhaseBuffer if it can acquire a TryLock. returns true if successful.
        bool tryToReadVisPhases(std::queue<double>& other);

        // for visualizer continuous channel
        void saveCustomChannelParametersToXml(XmlElement* channelElement, int channelNumber, InfoObjectCommon::InfoObjectType channelType) override;
        void loadCustomChannelParametersFromXml(XmlElement* channelElement, InfoObjectCommon::InfoObjectType channelType) override;

        /*
        * Circular distance between angles x and ref (in radians). The "cutoff" is the maximum
        * possible positive output; greater differences will be wrapped to negative angles.
        * For interpolation, and also RosePlot visualization.
        */
        static double circDist(double x, double ref, double cutoff = 2 * Dsp::doublePi);

    private:

        // ---- methods ----

        // Allow responding to stim events if a stimEventChannel is selected.
        void handleEvent(const EventChannel* eventInfo, const MidiMessage& event,
            int samplePosition = 0) override;

        // Sets frequency band (which resets lowCut and highCut to defaults for this band)
        void setBand(Band newBand, bool force = false);

        // Resets lowCut and highCut to defaults for the current band
        void resetCutsToDefaults();

        // Sets lowCut (which in turn influences highCut)
        void setLowCut(float newLowCut);

        // Sets highCut (which in turn influences lowCut)
        void setHighCut(float newHighCut);

        // Sets visContinuousChannel and updates the visualization filter
        void setVisContChan(int newChan);

        // Update the htScaleFactor depending on the lowCut and highCut of the filter.
        void updateScaleFactor();

        // Do glitch unwrapping
        void unwrapBuffer(float* wp, int nSamples, float lastPhase);

        // Do start-of-buffer smoothing
        void smoothBuffer(float* wp, int nSamples, float lastPhase);

        // Update subProcessorMap
        void updateSubProcessorMap();

        // Create an extra output channel for each processed input channel if PH_AND_MAG is selected
        void updateExtraChannels();

        // Deselect given channel in the "PARAMS" channel selector. Useful to ensure "extra channels"
        // remain deselected (so that they don't become active inputs if the # of inputs changes).
        // If "warn" is true, sends a warning that the channel was deselected due to incompatible
        // sample rate.
        void deselectChannel(int chan, bool warn);

        // Calls deselectChannel on each channel that is not currently an input. Only relevant
        // when the output mode is phase and magnitude.
        void deselectAllExtraChannels();

        /*
        * Check the visualization timestamp queue, clear any that are expired
        * (too late to calculate phase), and calculate phase of any that are ready.
        * sdbEndTs = timestamp 1 past end of current buffer
        * Precondition: chan is a valid input index.
        */
        void calcVisPhases(ActiveChannelInfo* acInfo, juce::int64 sdbEndTs);

        /*
        * Convenience method to call "update" on all channel info objects (which also updates
        * corresponding active channel info)
        */
        void updateAllChannels();

        /*
        * Convenience method to call "update" on all active channel info structs in the map.
        */
        void updateActiveChannels();

        /*
        * Enables an input channel to be processed. Returns false if this is prohibited due to
        * the channel's sample rate. If the channel is already active, does nothing.
        * Precondition: the channel passed in is a valid input channel index.
        */
        bool activateInputChannel(int chan);

        void deactivateInputChannel(int chan);

        // ---- static utility methods ----

        /*
        * arPredict: use autoregressive model of order to predict future data.
        *
        * lastSample points to the most recent sample of past data that will be used to
        * compute phase, and there must be at least stride * (order - 1) samples
        * preceding it in order to do the AR prediction.
        *
        * Input params is an array of coefficients of an AR model of length 'order'.
        *
        * Writes samps future data values to prediction.
        */
        static void arPredict(const ReverseStack& history, int interpCountdown, double* prediction,
            const double* params, int samps, int stride, int order);

        // Get the htScaleFactor for the given band's Hilbert transformer,
        // over the range from lowCut and highCut. This is the reciprocal of the geometric
        // mean (i.e. mean in decibels) of the maximum and minimum magnitude responses over the range.
        static double getScaleFactor(Band band, double lowCut, double highCut);

        // Execute the hilbert transformer on one sample and update the state.
        static double htFilterSamp(double input, Band band, Array<double>& state);

        // ---- customizable parameters ------

        // time to wait between AR model recalculations in ms
        int calcInterval;

        // order of the AR model
        int arOrder;

        OutputMode outputMode;

        // frequency band (determines which Hilbert transformer to use)
        Band band;

        // filter passband
        float highCut;
        float lowCut;

        // event channel to watch for phases to plot on the canvas (-1 = none)
        int visEventChannel;

        // channel to calculate phases from at received stim event times
        int visContinuousChannel;

        // ---- internals -------

        OwnedArray<ChannelInfo> channelInfo;

        // storage areas
        Array<double, CriticalSection> localARParams;
        Array<double> predSamps;
        Array<double> htTempState;
        Array<int> htInds;
        Array<std::complex<double>> htOutput;

        // approximate multiplier for the imaginary component output of the HT (depends on filter band)
        double htScaleFactor;

        // maps full source subprocessor IDs of incoming streams to indices of
        // corresponding subprocessors created here.
        HashMap<int, uint16> subProcessorMap;

        // delayed analysis for visualization

        // holds stimulation timestamps until the delayed phase is ready to be calculated
        std::queue<juce::int64> visTsBuffer;

        // for phases of stimulations, to be read by the visualizer.
        std::queue<double> visPhaseBuffer;
        CriticalSection visPhaseBufferCS;  // avoid race conditions when updating visualizer

        // event channel to send visualized phases over
        EventChannel* visPhaseChannel;

        JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(Node);
    };

}


#endif // PHASE_CALCULATOR_H_INCLUDED
