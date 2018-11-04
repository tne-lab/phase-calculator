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
#include <DspLib/Dsp.h>  // Filtering
#include <queue>
#include <array>

#include "FFTWWrapper.h"   // Fourier transform
#include "ARModeler.h"     // Autoregressive modeling
#include "HTransformers.h" // Hilbert transformers & frequency bands

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

/* each continuous channel has three possible states while acquisition is running:

    - NOT_FULL:     Not enough samples have arrived to fill the history fifo for this channel.
                    Wait for more  samples before calculating the autoregressive model parameters.

    - FULL_NO_AR:   The history fifo for this channel is now full, but AR parameters have not been calculated yet.
                    Tells the AR thread that it can start calculating the model and the main thread that it should still output zeros.

    - FULL_AR:      The history fifo is full and AR parameters have been calculated at least once.
                    In this state, the main thread uses the parameters to predict the future signal and output and calculate the phase.
*/
enum ChannelState { NOT_FULL, FULL_NO_AR, FULL_AR };

// Output mode - corresponds to itemIDs on the ComboBox
enum OutputMode { PH = 1, MAG, PH_AND_MAG, IM };

class PhaseCalculator : public GenericProcessor, public Thread
{
    friend class PhaseCalculatorEditor;

    // ------- constants ------------

    // default passband width if pushing lowCut down or highCut up to fix invalid range,
    // and also the minimum for lowCut.
    // (meant to be larger than actual minimum floating-point eps)
    static const float PASSBAND_EPS;

    // priority of the AR model calculating thread (0 = lowest, 10 = highest)
    static const int AR_PRIORITY = 3;

    // "glitch limit" (how long of a segment is allowed to be unwrapped or smoothed, in samples)
    static const int GLITCH_LIMIT = 200;

    // process length for real-time visualization ground truth hilbert transform
    // based on evaluation of phase error compared to offline processing
    static const int VIS_HILBERT_LENGTH = 65536;
    static const int VIS_TS_MIN_DELAY = 40000;
    static const int VIS_TS_MAX_DELAY = 48000;

public:

    PhaseCalculator();
    ~PhaseCalculator();

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
    Array<int> getActiveInputs();

    // ----- to create new channels for multiple outputs -------
    bool isGeneratesTimestamps() const override;
    int getNumSubProcessors() const override;
    float getSampleRate(int subProcessorIdx = 0) const override;
    float getBitVolts(int subProcessorIdx = 0) const override;

    // miscellaneous helper
    int getFullSourceId(int chan);

    // for the visualizer
    std::queue<double>& getVisPhaseBuffer(ScopedPointer<ScopedLock>& lock);

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

    // Sets arOrder (which in turn influences historyLength and the arModeler)
    void setAROrder(int newOrder);

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

    // Update historyLength to be the minimum possible size (depending on
    // VIS_HILBERT_LENGTH, hilbertLength, predictionLength, and arOrder).
    void updateHistoryLength();

    // Update the htScaleFactor depending on the lowCut and highCut of the filter.
    void updateScaleFactor();

    // Update the filters of active channels. From FilterNode code.
    void setFilterParameters();

    // If given input channel has an incompatible sample rate, deselect it (and send a warning).
    // Otherwise, save the downsampling multiplier and initialize downsampling phase.
    // Returns the new selection state of the channel.
    bool validateSampleRate(int chan);

    // Allocate memory for a new active input channel
    void addActiveChannel();

    // Do glitch unwrapping
    void unwrapBuffer(float* wp, int nSamples, int chan);

    // Do start-of-buffer smoothing
    void smoothBuffer(float* wp, int nSamples, int chan);

    // Update subProcessorMap
    void updateSubProcessorMap();

    // Create an extra output channel for each processed input channel if PH_AND_MAG is selected
    void updateExtraChannels();

    // Deselect given channel in the "PARAMS" channel selector. Useful to ensure "extra channels"
    // remain deselected (so that they don't become active inputs if the # of inputs changes).
    void deselectChannel(int chan);

    // Calls deselectChannel on each channel that is not currently an input. Only relevant
    // when the output mode is phase and magnitude.
    void deselectAllExtraChannels();

    /*
     * Check the visualization timestamp queue, clear any that are expired
     * (too late to calculate phase), and calculate phase of any that are ready.
     * sdbEndTs = timestamp 1 past end of current buffer
     */
    void calcVisPhases(juce::int64 sdbEndTs);

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
    static void arPredict(const double* lastSample, double* prediction,
        const double* params, int samps, int stride, int order);

    /*
     * hilbertManip: Hilbert transforms data in the frequency domain (including normalization by length of data).
     * Modifies fftData in place.
     */
    static void hilbertManip(FFTWArray* fftData);

    // Get the htScaleFactor for the given band's Hilbert transformer,
    // over the range from lowCut and highCut. This is the reciprocal of the geometric
    // mean (i.e. mean in decibels) of the maximum and minimum magnitude responses over the range.
    static double getScaleFactor(Band band, double lowCut, double highCut);

    // Execute the hilbert transformer on one sample and update the state.
    static double htFilterSamp(double input, Band band, Array<double>& state);

    // ---- customizable parameters ------

    // size of historyBuffer ( >= VIS_HILBERT_LENGTH and long enough to train AR model effectively)
    int historyLength;

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

    int numActiveChansAllocated = 0;

    // Storage area for filtered data to be used by the side thread to calculate AR model parameters
    // and by the visualization event handler to calculate acccurate phases at past event times.
    AudioBuffer<double> historyBuffer;

    // Keep track of how much of the historyBuffer is empty (per channel)
    Array<int> bufferFreeSpace;

    // Keeps track of each channel's state (see enum definition above)
    Array<ChannelState> chanState;

    // mutexes for sharedDataBuffer arrays, which are used in the side thread to calculate AR parameters.
    // since the side thread only READS sharedDataBuffer, only needs to be locked in the main thread when it's WRITTEN to.
    OwnedArray<CriticalSection> historyLock;

    // for autoregressive parameter calculation
    OwnedArray<ARModeler> arModelers;

    // for active input channels, equals the sample rate divided by HT_FS
    Array<int> sampleRateMultiple;

    // AR model parameters
    OwnedArray<Array<double>> arParams;
    OwnedArray<CriticalSection> arParamLock;

    // storage area for predicted samples (to compensate for group delay)
    Array<double> predSamps;

    // state of each hilbert transformer (persistent and temporary)
    OwnedArray<Array<double>> htState;
    Array<double> htTempState;

    // store indices of current buffer to feed into transformer
    Array<int> htInds;

    // store output of hilbert transformer
    Array<std::complex<double>> htOutput;

    // approximate multiplier for the imaginary component output of the HT (depends on filter band)
    double htScaleFactor;

    // keep track of each active channel's last non-interpolated ("computed") transformer output
    // and the number of samples by which this sample precedes the start of the next buffer
    // (ranges from 1 to sampleRateMultiple).
    Array<std::complex<double>> lastComputedSample;
    Array<int> dsOffset;

    // keep track of last phase output, for glitch correction
    Array<float> lastPhase;

    // maps full source subprocessor IDs of incoming streams to indices of
    // corresponding subprocessors created here.
    HashMap<int, uint16> subProcessorMap;

    // delayed analysis for visualization
    FFTWArray visHilbertBuffer;
    FFTWPlan  visForwardPlan, visBackwardPlan;

    // holds stimulation timestamps until the delayed phase is ready to be calculated
    std::queue<juce::int64> visTsBuffer;

    // for phases of stimulations, to be read by the visualizer.
    std::queue<double> visPhaseBuffer;
    CriticalSection visPhaseBufferLock;  // avoid race conditions when updating visualizer

    // event channel to send visualized phases over
    EventChannel* visPhaseChannel;

    // filter design copied from FilterNode
    typedef Dsp::SmoothedFilterDesign
        <Dsp::Butterworth::Design::BandPass // design type
        <2>,                                // order
        1,                                  // number of channels
        Dsp::DirectFormII>                  // realization
        BandpassFilterBase;

    class BandpassFilter : public BandpassFilterBase
    {
    public:
        BandpassFilter() : BandpassFilterBase(1) {} // # of transition samples
    };

    OwnedArray<BandpassFilter> filters;
    BandpassFilter visReverseFilter;

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(PhaseCalculator);
};

#endif // PHASE_CALCULATOR_H_INCLUDED
