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
#include "FFTWWrapper.h" // Fourier transform
#include "ARModeler.h"   // Autoregressive modeling

// parameter indices
enum Param
{
    PRED_LENGTH,
    RECALC_INTERVAL,
    AR_ORDER,
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
    friend class ProcessBufferSlider;
public:

    PhaseCalculator();
    ~PhaseCalculator();

    bool hasEditor() const override;

    AudioProcessorEditor* createEditor() override;

    void setParameter(int parameterIndex, float newValue) override;

    void process(AudioSampleBuffer& buffer) override;

    bool enable() override;
    bool disable() override;

    // calculate the fraction of the processing length is AR-predicted data points
    // (i.e. predictionLength / hilbertLength, as a float)
    float getPredictionRatio();

    // thread code - recalculates AR parameters.
    void run() override;

    // handle changing number of channels
    void updateSettings() override;

    // ----- to create new channels for multiple outputs -------
    bool isGeneratesTimestamps() const override;
    int getNumSubProcessors() const override;
    float getSampleRate(int subProcessorIdx = 0) const override;
    float getBitVolts(int subProcessorIdx = 0) const override;

    // for the visualizer
    std::queue<double>& getVisPhaseBuffer(ScopedPointer<ScopedLock>& lock);

    // for visualizer continuous channel
    void saveCustomChannelParametersToXml(XmlElement* channelElement, int channelNumber, InfoObjectCommon::InfoObjectType channelType) override;
    void loadCustomChannelParametersFromXml(XmlElement* channelElement, InfoObjectCommon::InfoObjectType channelType) override;

private:

    // ---- methods ----

    // Allow responding to stim events if a stimEventChannel is selected.
    void handleEvent(const EventChannel* eventInfo, const MidiMessage& event,
        int samplePosition = 0) override;

    // Set hilbertLength and predictionLength, reallocating fields that depend on these.
    void setHilbertAndPredLength(int newHilbertLength, int newPredictionLength);

    // Update historyLength to be the minimum possible size (depending on
    // VIS_HILBERT_LENGTH, hilbertLength, predictionLength, and arOrder).
    void updateHistoryLength();

    // Update minimum Nyquist frequency of active channels based on sample rates
    void updateMinNyquist();

    // Update the filters of active channels. From FilterNode code.
    void setFilterParameters();

    // Do glitch unwrapping
    void unwrapBuffer(float* wp, int nSamples, int chan);

    // Do start-of-buffer smoothing
    void smoothBuffer(float* wp, int nSamples, int chan);

    // Update subProcessorMap
    void updateSubProcessorMap();

    // Create an extra output channel for each processed input channel if PH_AND_MAG is selected
    void updateExtraChannels();

    /* 
     * Check the visualization timestamp queue, clear any that are expired
     * (too late to calculate phase), and calculate phase of any that are ready.
     * sdbEndTs = timestamp 1 past end of current buffer
     */
    void calcVisPhases(juce::int64 sdbEndTs);

    // ---- static utility methods ----

    /*
     * arPredict: use autoregressive model of order to predict future data.
     * Input params is an array of coefficients of an AR model of length 'order'.
     * Writes writeNum future data values starting at location writeStart.
     * *** assumes there are at least 'order' existing data points *before* writeStart
     * to use to calculate future data points.
     */
    static void arPredict(double* writeStart, int writeNum, const double* params, int order);

    /*
     * hilbertManip: Hilbert transforms data in the frequency domain (including normalization by length of data).
     * Modifies fftData in place.
     */
    static void hilbertManip(FFTWArray* fftData);

    // ---- customizable parameters ------
    
    // number of samples to pass through the Hilbert transform in the main calculation
    int hilbertLength;

    // number of future values to predict (0 <= predictionLength <= hilbertLength - AR_ORDER)
    int predictionLength;

    // size of historyBuffer ( = max(VIS_HILBERT_LENGTH, hilbertLength - predictionLength))
    int historyLength;

    // time to wait between AR model recalculations in ms
    int calcInterval;

    // order of the AR model
    int arOrder;

    OutputMode outputMode;

    // filter passband
    double highCut;
    double lowCut;

    // event channel to watch for phases to plot on the canvas (-1 = none)
    int visEventChannel;

    // channel to calculate phases from at received stim event times
    int visContinuousChannel;

    // ---- internals -------

    // Storage area for filtered data to be read by the main thread to fill hilbertBuffer,
    // by the side thread to calculate AR model parameters, and by the visualization event
    // handler to calculate acccurate phases at past event times.
    AudioBuffer<double> historyBuffer;

    // Keep track of how much of the historyBuffer is empty (per channel)
    Array<int> bufferFreeSpace;

    // Keeps track of each channel's state (see enum definition above)
    Array<ChannelState> chanState;

    // Buffers for FFTW processing
    OwnedArray<FFTWArray> hilbertBuffer;
    
    // Plans for the FFTW Fourier Transform library
    OwnedArray<FFTWPlan> forwardPlan;  // FFT
    OwnedArray<FFTWPlan> backwardPlan; // IFFT

    // mutexes for sharedDataBuffer arrays, which are used in the side thread to calculate AR parameters.
    // since the side thread only READS sharedDataBuffer, only needs to be locked in the main thread when it's WRITTEN to.
    OwnedArray<CriticalSection> historyLock;

    // for autoregressive parameter calculation
    ARModeler arModeler;

    // keeps track of each channel's last output sample, to be used for smoothing.
    Array<float> lastSample;

    // AR model parameters
    OwnedArray<Array<double>> arParams;
    OwnedArray<CriticalSection> arParamLock;

    // so that the warning message only gets sent once per run
    bool haveSentWarning;

    // maps full source subprocessor IDs of incoming streams to indices of
    // corresponding subprocessors created here.
    HashMap<int, uint16> subProcessorMap;

    // for filtering - min Nyquist frequency over active channels
    // (upper limit for highCut)
    float minNyquist;

    // delayed analysis for visualization
    FFTWArray visHilbertBuffer;
    FFTWPlan  visForwardPlan, visBackwardPlan;

    // holds stimulation timestamps until the delayed phase is ready to be calculated
    std::queue<juce::int64> visTsBuffer;

    // for phases of stimulations, to be read by the visualizer.
    std::queue<double> visPhaseBuffer;
    CriticalSection visPhaseBufferLock;  // avoid race conditions when updating visualizer

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

    // -------static------------

    // priority of the AR model calculating thread (0 = lowest, 10 = highest)
    static const int AR_PRIORITY = 3;

    // "glitch limit" (how long of a segment is allowed to be unwrapped or smoothed, in samples)
    static const int GLITCH_LIMIT = 200;

    // process length limits (powers of 2)
    static const int MIN_PLEN_POW = 9;
    static const int MAX_PLEN_POW = 16;

    // process length for real-time visualization ground truth hilbert transform
    // based on evaluation of phase error compared to offline processing
    static const int VIS_HILBERT_LENGTH = 65536;
    static const int VIS_TS_MIN_DELAY = 40000;
    static const int VIS_TS_MAX_DELAY = 48000;

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(PhaseCalculator);
};

// timer class for use when calculating AR parameters
class ARTimer : public Timer
{
public:
    ARTimer();
    ~ARTimer();
    void timerCallback();
    // Returns whether hasRung is true and resets it to false.
    bool check();
private:
    // True if timer has reached 0 since last time it was checked.
    bool hasRung;
};

#endif // PHASE_CALCULATOR_H_INCLUDED