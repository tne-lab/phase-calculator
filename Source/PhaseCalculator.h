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
#include <DspLib.h>           // Filtering
#include <OpenEphysFFTW.h>    // Fourier transform

#include <queue>
#include <utility>  // pair

#include "ARModeler.h"     // Autoregressive modeling
#include "HTransformers.h" // Hilbert transformers & frequency bands

namespace PhaseCalculator
{
    // forward declarations
    struct ChannelInfo;
    class Node;

    
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
        ActiveChannelInfo(const ChannelInfo* cInfo);

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

        const ChannelInfo* chanInfo;

        void waitForThreadToExit()
        {
            if (bufferResizeThread->isThreadRunning())
            {
                bufferResizeThread->waitForThreadToExit(1000);
            }
        }

    private:

        class BufferResizeThread : public Thread
        {
        public:
            BufferResizeThread(FFTWTransformableArray* buffer_) : Thread("buffer resize"),
                buffer(buffer_), bufferSize(0)
            {

            }

            ~BufferResizeThread() { 
            
                if (isThreadRunning())
                {
                    waitForThreadToExit(5000);
                }
            }

            void setSize(int lengthMultiplier)
            {
                bufferSize = lengthMultiplier;
            }

            void run()
            {
                buffer->resize(bufferSize);
            }

        private:
            FFTWTransformableArray* buffer;

            int bufferSize;
        };

        std::unique_ptr<BufferResizeThread> bufferResizeThread;

        JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(ActiveChannelInfo);
    };

    struct ChannelInfo
    {
        ChannelInfo(const DataStream* ds, int i);

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
        std::unique_ptr<ActiveChannelInfo> acInfo;
        const DataStream* stream;

    private:

        bool isActivated; 

        JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(ChannelInfo);
    };


    class Settings
    {
    public:
        Settings();

        ~Settings() { }

        // Returns array of active channels that only includes inputs (not extra outputs)
        Array<int> getActiveInputs() const;

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

        // Sets frequency band (which resets lowCut and highCut to defaults for this band)
        void setBand(Band newBand, bool force = false);

        // Resets lowCut and highCut to defaults for the current band
        void resetCutsToDefaults();

        // Update the htScaleFactor depending on the lowCut and highCut of the filter.
        void updateScaleFactor();

         // Get the htScaleFactor for the given band's Hilbert transformer,
        // over the range from lowCut and highCut. This is the reciprocal of the geometric
        // mean (i.e. mean in decibels) of the maximum and minimum magnitude responses over the range.
        static double getScaleFactor(Band band, double lowCut, double highCut);

        // approximate multiplier for the imaginary component output of the HT (depends on filter band)
        double htScaleFactor;

        OwnedArray<ChannelInfo> channelInfo;

        // ---- customizable parameters ------

        // time to wait between AR model recalculations in ms
        int calcInterval;

        // order of the AR model
        int arOrder;

        // frequency band (determines which Hilbert transformer to use)
        Band band;

        // filter passband
        float highCut;
        float lowCut;

        // event channel to watch for phases to plot on the canvas (-1 = none)
        int visEventChannel;

        // channel to calculate phases from at received stim event times
        int visContinuousChannel;

        Array<double> predSamps;
        Array<double> htTempState;
    };

    class Node : public GenericProcessor, public Thread
    {
        friend class Editor;
    public:

        /** Constructor */
        Node();

        /** Destructor */
        ~Node() { }

        /** Creates the PhaseCalculatorEditor */
        AudioProcessorEditor* createEditor() override;

        /** 
            Overwrites a subset of channels with the instantaneous phase in a given
            frequency band.
         */
        void process(AudioBuffer<float>& buffer) override;

        /** Starts phase calculation thread */
        bool startAcquisition() override;

        /** Stops phase calculation code */
        bool stopAcquisition() override;

        /** thread code - recalculates AR parameters. */
        void run() override;

        /** Called whenever inputs are changed  */
        void updateSettings() override;

        /** Called whenever a parameter's value is changed (called by GenericProcessor::setParameter())*/
        void parameterValueChanged(Parameter* param) override;

        /** reads from the visPhaseBuffer if it can acquire a TryLock. returns true if successful. */
        bool tryToReadVisPhases(std::queue<double>& other);

        /** Returns array of active channels that only includes inputs (not extra outputs) */
        Array<int> getActiveChannels();

        /*
        * Circular distance between angles x and ref (in radians). The "cutoff" is the maximum
        * possible positive output; greater differences will be wrapped to negative angles.
        * For interpolation, and also RosePlot visualization.
        */
        static double circDist(double x, double ref, double cutoff = 2 * Dsp::doublePi);

        /** Get the current selected stream */
        juce::uint16 getSelectedStream() { return selectedStream; }

        /** Set the current selected stream */
        void setSelectedStream(juce::uint16 streamId);

    private:

        // ---- methods ----

        /** Responds to incoming events if a stimEventChannel is selected. */
        void handleTTLEvent(TTLEventPtr event) override;

        /** Perform glitch unwrapping */
        void unwrapBuffer(float* wp, int nSamples, float lastPhase);

        /** Do start-of-buffer smoothing */
        void smoothBuffer(float* wp, int nSamples, float lastPhase);

        /*
        * Check the visualization timestamp queue, clear any that are expired
        * (too late to calculate phase), and calculate phase of any that are ready.
        * sdbEndTs = timestamp 1 past end of current buffer
        * Precondition: chan is a valid input index.
        */
        void calcVisPhases(ActiveChannelInfo* acInfo, juce::int64 sdbEndTs);

        /** Sets visContinuousChannel and updates the visualization filter */
        void setVisContChan(int newChan);
        
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


        /** Execute the hilbert transformer on one sample and update the state. */
        static double htFilterSamp(double input, Band band, Array<double>& state);

        // ---- internals -------

        StreamSettings<Settings> settings;

        /** Only one stream can be activated at a time*/
        uint16 selectedStream;

        // storage areas
        Array<double, CriticalSection> localARParams;
        Array<int> htInds;
        Array<std::complex<double>> htOutput;

        // delayed analysis for visualization

        // holds stimulation timestamps until the delayed phase is ready to be calculated
        std::queue<juce::int64> visTsBuffer;

        // for phases of stimulations, to be read by the visualizer.
        std::queue<double> visPhaseBuffer;
        CriticalSection visPhaseBufferCS;  // avoid race conditions when updating visualizer

        /** Notify Node thread to update it's list of active channels and find maximum history length */
        bool activeChansNeedsUpdate;


        JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(Node);
    };

}


#endif // PHASE_CALCULATOR_H_INCLUDED
