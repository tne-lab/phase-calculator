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


#ifndef PHASE_CALCULATOR_EDITOR_H_INCLUDED
#define PHASE_CALCULATOR_EDITOR_H_INCLUDED

#include <VisualizerEditorHeaders.h>

#include "PhaseCalculator.h"

class PhaseCalculatorEditor
    : public VisualizerEditor
    , public ComboBox::Listener
    , public Label::Listener
{
    friend class RosePlot;  // to access label updating method
public:
    PhaseCalculatorEditor(PhaseCalculator* parentNode, bool useDefaultParameterEditors = false);
    ~PhaseCalculatorEditor();

    // implements ComboBox::Listener
    void comboBoxChanged(ComboBox* comboBoxThatHasChanged) override;

    // implements Label::Listener
    void labelTextChanged(Label* labelThatHasChanged) override;

    // overrides GenericEditor
    void sliderEvent(Slider* slider) override;

    // update display based on current channel
    void channelChanged(int chan, bool newState) override;

    // enable/disable controls when acquisiton starts/ends
    void startAcquisition() override;
    void stopAcquisition() override;

    Visualizer* createNewCanvas() override;

    void updateSettings() override;

    void saveCustomParameters(XmlElement* xml) override;
    void loadCustomParameters(XmlElement* xml) override;

    // display updaters - do not trigger listeners.
    void refreshLowCut();
    void refreshHighCut();
    void refreshPredLength();
    void refreshHilbertLength();
    void refreshVisContinuousChan();

private:

    // Utilities for parsing entered values

    template<typename T>
    static T fromString(const char* in);

    /* 
     * Return whether the control contained a valid input between min and max, inclusive.
     * If so, it is stored in *out and the control is updated with the parsed input.
     * Otherwise, the control is reset to defaultValue.
     *
     * In header to make sure specializations not used in PhaseCalculatorEditor.cpp
     * are still available to other translation units.
     */

    template<typename Ctrl, typename T>
    static bool updateControl(Ctrl* c, const T min, const T max,
        const T defaultValue, T* out)
    {
        T parsedVal;
        try
        {
            parsedVal = fromString<T>(c->getText().toRawUTF8());
        }
        catch (const std::logic_error&)
        {
            c->setText(String(defaultValue), dontSendNotification);
            return false;
        }

        *out = jmax(min, jmin(max, parsedVal));

        c->setText(String(*out), dontSendNotification);
        return true;
    }

    // keep track of the record status of each "extra" channel
    // this is all a bit kludgy, but also temporary until creating continuous
    // channels in non-source processors is officially implemented.
    class ExtraChanManager : public Button::Listener
    {
    public:
        ExtraChanManager(const PhaseCalculator* processor);

        // keeps recordStatus in sync with extra channel record buttons.
        void buttonClicked(Button* button) override;

        // adds or removes recordStatus entry for extra chan
        // corresponding to given input chan.
        void addExtraChan(int inputChan, const Array<int>& activeInputs);
        void removeExtraChan(int inputChan, const Array<int>& activeInputs);
        void resize(int numExtraChans);

        bool getRecordStatus(int extraChan) const;
    private:
        const PhaseCalculator* p;
        Array<bool> recordStatus;
    };

    int prevExtraChans;
    ExtraChanManager extraChanManager;

    ScopedPointer<Label>    lowCutLabel;
    ScopedPointer<Label>    lowCutEditable;
    ScopedPointer<Label>    highCutLabel;
    ScopedPointer<Label>    highCutEditable;
    
    ScopedPointer<Label>    hilbertLengthLabel;
    ScopedPointer<Label>    hilbertLengthUnitLabel;
    ScopedPointer<ComboBox> hilbertLengthBox;

    LookAndFeel_V3        v3LookAndFeel;
    ScopedPointer<Slider> predLengthSlider;
    ScopedPointer<Label>  pastLengthLabel;
    ScopedPointer<Label>  pastLengthEditable;
    ScopedPointer<Label>  predLengthLabel;
    ScopedPointer<Label>  predLengthEditable;

    ScopedPointer<Label>    recalcIntervalLabel;
    ScopedPointer<Label>    recalcIntervalEditable;
    ScopedPointer<Label>    recalcIntervalUnit;

    ScopedPointer<Label>    arOrderLabel;
    ScopedPointer<Label>    arOrderEditable;

    ScopedPointer<Label>    outputModeLabel;
    ScopedPointer<ComboBox> outputModeBox;

    // constants

    const String HILB_LENGTH_TOOLTIP = "Change the total amount of data used to calculate the phase (powers of 2 are best)";
    const String PRED_LENGTH_TOOLTIP = "Select how much past vs. predicted data to use when calculating the phase";
    const String RECALC_INTERVAL_TOOLTIP = "Time to wait between calls to update the autoregressive models";
    const String AR_ORDER_TOOLTIP = "Order of the autoregressive models used to predict future data";
    const String OUTPUT_MODE_TOOLTIP = "Which component of the analytic signal to output. If 'PH+MAG' is selected, " +
        String("creates a second channel for each enabled input channel and outputs phases ") +
        "on the original channels and magnitudes on the corresponding new channels.";

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(PhaseCalculatorEditor);
};

#endif // PHASE_CALCULATOR_EDITOR_H_INCLUDED