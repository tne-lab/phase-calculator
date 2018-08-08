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

#include "PhaseCalculatorEditor.h"
#include "PhaseCalculatorCanvas.h"
#include <climits> // INT_MAX

const float PhaseCalculatorEditor::PASSBAND_EPS = 0.01F;

PhaseCalculatorEditor::PhaseCalculatorEditor(GenericProcessor* parentNode, bool useDefaultParameterEditors)
    : VisualizerEditor     (parentNode, 325, useDefaultParameterEditors)
{
    tabText = "Event Phase Plot";
    int filterWidth = 80;

    PhaseCalculator* processor = static_cast<PhaseCalculator*>(parentNode);

    lowCutLabel = new Label("lowCutL", "Low cut");
    lowCutLabel->setBounds(10, 30, 80, 20);
    lowCutLabel->setFont(Font("Small Text", 12, Font::plain));
    lowCutLabel->setColour(Label::textColourId, Colours::darkgrey);
    addAndMakeVisible(lowCutLabel);

    lowCutEditable = new Label("lowCutE");
    lowCutEditable->setEditable(true);
    lowCutEditable->addListener(this);
    lowCutEditable->setBounds(15, 47, 60, 18);
    lowCutEditable->setText(String(processor->lowCut), dontSendNotification);
    lowCutEditable->setColour(Label::backgroundColourId, Colours::grey);
    lowCutEditable->setColour(Label::textColourId, Colours::white);
    addAndMakeVisible(lowCutEditable);

    highCutLabel = new Label("highCutL", "High cut");
    highCutLabel->setBounds(10, 70, 80, 20);
    highCutLabel->setFont(Font("Small Text", 12, Font::plain));
    highCutLabel->setColour(Label::textColourId, Colours::darkgrey);
    addAndMakeVisible(highCutLabel);

    highCutEditable = new Label("highCutE");
    highCutEditable->setEditable(true);
    highCutEditable->addListener(this);
    highCutEditable->setBounds(15, 87, 60, 18);
    highCutEditable->setText(String(processor->highCut), dontSendNotification);
    highCutEditable->setColour(Label::backgroundColourId, Colours::grey);
    highCutEditable->setColour(Label::textColourId, Colours::white);
    addAndMakeVisible(highCutEditable);

    hilbertLengthLabel = new Label("hilbertLength", "Buffer length:");
    hilbertLengthLabel->setBounds(filterWidth + 8, 25, 180, 20);
    hilbertLengthLabel->setFont(Font("Small Text", 12, Font::plain));
    hilbertLengthLabel->setColour(Label::textColourId, Colours::darkgrey);
    addAndMakeVisible(hilbertLengthLabel);

    hilbertLengthBox = new ComboBox("Buffer size");
    hilbertLengthBox->setEditableText(true);
    for (int pow = PhaseCalculator::MIN_PLEN_POW; pow <= PhaseCalculator::MAX_PLEN_POW; ++pow)
    {
        hilbertLengthBox->addItem(String(1 << pow), pow);
    }
    hilbertLengthBox->setText(String(processor->hilbertLength), dontSendNotification);
    hilbertLengthBox->setTooltip(HILB_LENGTH_TOOLTIP);
    hilbertLengthBox->setBounds(filterWidth + 10, 45, 80, 20);
    hilbertLengthBox->addListener(this);
    addAndMakeVisible(hilbertLengthBox);

    hilbertLengthUnitLabel = new Label("hilbertLengthUnit", "Samp.");
    hilbertLengthUnitLabel->setBounds(filterWidth + 90, 45, 40, 20);
    hilbertLengthUnitLabel->setFont(Font("Small Text", 12, Font::plain));
    hilbertLengthUnitLabel->setColour(Label::textColourId, Colours::darkgrey);
    addAndMakeVisible(hilbertLengthUnitLabel);

    pastLengthLabel = new Label("pastLengthL", "Past:");
    pastLengthLabel->setBounds(filterWidth + 8, 85, 60, 15);
    pastLengthLabel->setFont(Font("Small Text", 12, Font::plain));
    pastLengthLabel->setColour(Label::backgroundColourId, Colour(230, 168, 0));
    pastLengthLabel->setColour(Label::textColourId, Colours::darkgrey);
    addAndMakeVisible(pastLengthLabel);

    predLengthLabel = new Label("predLengthL", "Future:");
    predLengthLabel->setBounds(filterWidth + 70, 85, 60, 15);
    predLengthLabel->setFont(Font("Small Text", 12, Font::plain));
    predLengthLabel->setColour(Label::backgroundColourId, Colour(102, 140, 255));
    predLengthLabel->setColour(Label::textColourId, Colours::darkgrey);
    addAndMakeVisible(predLengthLabel);

    pastLengthEditable = new Label("pastLengthE");
    pastLengthEditable->setEditable(true);
    pastLengthEditable->addListener(this);
    pastLengthEditable->setBounds(filterWidth + 8, 102, 60, 18);
    pastLengthEditable->setColour(Label::backgroundColourId, Colours::grey);
    pastLengthEditable->setColour(Label::textColourId, Colours::white);

    predLengthEditable = new Label("predLengthE");
    predLengthEditable->setEditable(true);
    predLengthEditable->addListener(this);
    predLengthEditable->setBounds(filterWidth + 70, 102, 60, 18);
    predLengthEditable->setColour(Label::backgroundColourId, Colours::grey);
    predLengthEditable->setColour(Label::textColourId, Colours::white);

    predLengthSlider = new ProcessBufferSlider("predLength");
    predLengthSlider->setBounds(filterWidth + 8, 70, 122, 10);
    predLengthSlider->setColour(Slider::thumbColourId, Colour(255, 187, 0));
    predLengthSlider->setColour(Slider::backgroundColourId, Colour(51, 102, 255));
    predLengthSlider->setTooltip(PRED_LENGTH_TOOLTIP);
    predLengthSlider->addListener(this);
    predLengthSlider->updateFromProcessor(processor);
    addAndMakeVisible(predLengthSlider);
    addAndMakeVisible(pastLengthEditable);
    addAndMakeVisible(predLengthEditable);

    recalcIntervalLabel = new Label("recalcL", "AR Refresh:");
    recalcIntervalLabel->setBounds(filterWidth + 140, 25, 100, 20);
    recalcIntervalLabel->setFont(Font("Small Text", 12, Font::plain));
    recalcIntervalLabel->setColour(Label::textColourId, Colours::darkgrey);
    addAndMakeVisible(recalcIntervalLabel);

    recalcIntervalEditable = new Label("recalcE");
    recalcIntervalEditable->setEditable(true);
    recalcIntervalEditable->addListener(this);
    recalcIntervalEditable->setBounds(filterWidth + 145, 44, 55, 18);
    recalcIntervalEditable->setColour(Label::backgroundColourId, Colours::grey);
    recalcIntervalEditable->setColour(Label::textColourId, Colours::white);
    recalcIntervalEditable->setText(String(processor->calcInterval), dontSendNotification);
    recalcIntervalEditable->setTooltip(RECALC_INTERVAL_TOOLTIP);
    addAndMakeVisible(recalcIntervalEditable);

    recalcIntervalUnit = new Label("recalcU", "ms");
    recalcIntervalUnit->setBounds(filterWidth + 200, 47, 25, 15);
    recalcIntervalUnit->setFont(Font("Small Text", 12, Font::plain));
    recalcIntervalUnit->setColour(Label::textColourId, Colours::darkgrey);
    addAndMakeVisible(recalcIntervalUnit);

    arOrderLabel = new Label("arOrderL", "Order:");
    arOrderLabel->setBounds(filterWidth + 140, 65, 60, 20);
    arOrderLabel->setFont(Font("Small Text", 12, Font::plain));
    arOrderLabel->setColour(Label::textColourId, Colours::darkgrey);
    addAndMakeVisible(arOrderLabel);

    arOrderEditable = new Label("arOrderE");
    arOrderEditable->setEditable(true);
    arOrderEditable->addListener(this);
    arOrderEditable->setBounds(filterWidth + 195, 66, 25, 18);
    arOrderEditable->setColour(Label::backgroundColourId, Colours::grey);
    arOrderEditable->setColour(Label::textColourId, Colours::white);
    arOrderEditable->setText(String(processor->arOrder), sendNotificationAsync);
    arOrderEditable->setTooltip(AR_ORDER_TOOLTIP);
    addAndMakeVisible(arOrderEditable);

    outputModeLabel = new Label("outputModeL", "Output:");
    outputModeLabel->setBounds(filterWidth + 140, 87, 70, 20);
    outputModeLabel->setFont(Font("Small Text", 12, Font::plain));
    outputModeLabel->setColour(Label::textColourId, Colours::darkgrey);
    addAndMakeVisible(outputModeLabel);

    outputModeBox = new ComboBox("outputModeB");
    outputModeBox->addItem("PHASE", PH);
    outputModeBox->addItem("MAG", MAG);
    outputModeBox->addItem("PH+MAG", PH_AND_MAG);
    outputModeBox->addItem("IMAG", IM);
    outputModeBox->setSelectedId(processor->outputMode);
    outputModeBox->setTooltip(OUTPUT_MODE_TOOLTIP);
    outputModeBox->setBounds(filterWidth + 145, 105, 76, 19);
    outputModeBox->addListener(this);
    addAndMakeVisible(outputModeBox);

    // new channels should be disabled by default
    channelSelector->paramButtonsToggledByDefault(false);
}

PhaseCalculatorEditor::~PhaseCalculatorEditor() {}

void PhaseCalculatorEditor::comboBoxChanged(ComboBox* comboBoxThatHasChanged)
{
    PhaseCalculator* processor = static_cast<PhaseCalculator*>(getProcessor());

    if (comboBoxThatHasChanged == hilbertLengthBox)
    {
        int newId = hilbertLengthBox->getSelectedId();
        int newHilbertLength;
        if (newId) // one of the items in the list is selected
        {            
            newHilbertLength = (1 << newId);
        }
        else
        {
            // try to parse input
            String input = hilbertLengthBox->getText();
            int parsedInt;
            try
            {
                parsedInt = std::stoi(input.toRawUTF8());
            }
            catch (const std::logic_error&)
            {
                hilbertLengthBox->setText(String(processor->hilbertLength), dontSendNotification);
                return;
            }

            newHilbertLength = jmax(1 << PhaseCalculator::MIN_PLEN_POW,
                jmin(1 << PhaseCalculator::MAX_PLEN_POW, parsedInt));
           
            hilbertLengthBox->setText(String(newHilbertLength), dontSendNotification);
        }

        // update AR order if necessary
        if (newHilbertLength < processor->arOrder)
        {
            arOrderEditable->setText(String(newHilbertLength), sendNotificationSync);
            CoreServices::sendStatusMessage("AR order snapped to maximum value for this processing buffer length");
        }

        // calculate predLength
        float currRatio = processor->getPredictionRatio();
        float newPredLength;
        newPredLength = jmin(static_cast<int>(roundf(currRatio * newHilbertLength)), newHilbertLength - processor->arOrder);

        // change both at once
        processor->setHilbertAndPredLength(newHilbertLength, newPredLength);

        // update slider
        predLengthSlider->updateFromProcessor(processor);
    }
    else if (comboBoxThatHasChanged == outputModeBox)
    {
        processor->setParameter(OUTPUT_MODE, static_cast<float>(outputModeBox->getSelectedId()));
    }
}

void PhaseCalculatorEditor::labelTextChanged(Label* labelThatHasChanged)
{
    PhaseCalculator* processor = static_cast<PhaseCalculator*>(getProcessor());

    int sliderMin = static_cast<int>(predLengthSlider->getRealMinValue());
    int sliderMax = static_cast<int>(predLengthSlider->getMaximum());

    if (labelThatHasChanged == pastLengthEditable)
    {
        int intInput;
        bool valid = updateIntLabel(labelThatHasChanged, sliderMin, sliderMax,
            static_cast<int>(predLengthSlider->getValue()), &intInput);

        if (valid)
        {
            int newPredLength = sliderMax - intInput;
            predLengthSlider->setValue(intInput, dontSendNotification);
            predLengthEditable->setText(String(newPredLength), dontSendNotification);
            processor->setParameter(PRED_LENGTH, static_cast<float>(newPredLength));
        }        
    }
    else if (labelThatHasChanged == predLengthEditable)
    {
        int intInput;
        bool valid = updateIntLabel(labelThatHasChanged, 0, sliderMax - sliderMin,
            sliderMax - static_cast<int>(predLengthSlider->getValue()), &intInput);
        
        if (valid)
        {
            int newPastLength = sliderMax - intInput;
            predLengthSlider->setValue(newPastLength, dontSendNotification);
            pastLengthEditable->setText(String(newPastLength), dontSendNotification);
            processor->setParameter(PRED_LENGTH, static_cast<float>(intInput));
        }
    }
    else if (labelThatHasChanged == recalcIntervalEditable)
    {
        int intInput;
        bool valid = updateIntLabel(labelThatHasChanged, 0, INT_MAX, processor->calcInterval, &intInput);

        if (valid)
        {
            processor->setParameter(RECALC_INTERVAL, static_cast<float>(intInput));
        }
    }
    else if (labelThatHasChanged == arOrderEditable)
    {
        int intInput;
        bool valid = updateIntLabel(labelThatHasChanged, 1, processor->hilbertLength, processor->arOrder, &intInput);

        if (valid)
        {
            processor->setParameter(AR_ORDER, static_cast<float>(intInput));
            // update slider's minimum value, and predictionLength and historyLength if necessary
            predLengthSlider->updateFromProcessor(processor);
        }
    }
    else if (labelThatHasChanged == lowCutEditable)
    {
        float floatInput;
        bool valid = updateFloatLabel(labelThatHasChanged, PASSBAND_EPS, processor->minNyquist - PASSBAND_EPS,
            processor->lowCut, &floatInput);

        if (!valid) { return; }
        
        if (floatInput > processor->highCut)
        {   
            // push highCut up
            highCutEditable->setText(String(floatInput + PASSBAND_EPS), sendNotificationSync);
        }
        processor->setParameter(LOWCUT, floatInput);
    }
    else if (labelThatHasChanged == highCutEditable)
    {
        float floatInput;
        bool valid = updateFloatLabel(labelThatHasChanged, 2 * PASSBAND_EPS, processor->minNyquist,
            processor->highCut, &floatInput);

        if (!valid) { return; }

        if (floatInput < processor->lowCut)
        {
            // push lowCut down
            lowCutEditable->setText(String(floatInput - PASSBAND_EPS), sendNotificationSync);
        }
        processor->setParameter(HIGHCUT, floatInput);
    }
}

void PhaseCalculatorEditor::sliderEvent(Slider* slider)
{
    if (slider == predLengthSlider)
    {
        // At this point, a snapValue call has ensured that the new value is valid.
        int newVal = static_cast<int>(slider->getValue());
        int maxVal = static_cast<int>(slider->getMaximum());
        pastLengthEditable->setText(String(newVal), dontSendNotification);
        predLengthEditable->setText(String(maxVal - newVal), dontSendNotification);
        
        getProcessor()->setParameter(PRED_LENGTH, static_cast<float>(maxVal - newVal));
    }
}

void PhaseCalculatorEditor::channelChanged(int chan, bool newState)
{
    // If not an output channel, update signal chain. This should take care of:
    //     - adding/removing output channels if necessary
    //     - updating the minimum Nyquist frequency of active channels (and highCut if necessary)
    //     - updating the available continuous channels for visualizer
    if (chan < getProcessor()->getNumInputs())
    {
        CoreServices::updateSignalChain(this);
    }
}

void PhaseCalculatorEditor::startAcquisition()
{
    GenericEditor::startAcquisition();
    hilbertLengthBox->setEnabled(false);
    predLengthSlider->setEnabled(false);
    pastLengthEditable->setEnabled(false);
    predLengthEditable->setEnabled(false);
    lowCutEditable->setEnabled(false);
    highCutEditable->setEnabled(false);
    arOrderEditable->setEnabled(false);
    outputModeBox->setEnabled(false);
    channelSelector->inactivateButtons();
}

void PhaseCalculatorEditor::stopAcquisition()
{
    GenericEditor::stopAcquisition();
    hilbertLengthBox->setEnabled(true);
    predLengthSlider->setEnabled(true);
    pastLengthEditable->setEnabled(true);
    predLengthEditable->setEnabled(true);
    lowCutEditable->setEnabled(true);
    highCutEditable->setEnabled(true);
    arOrderEditable->setEnabled(true);
    outputModeBox->setEnabled(true);
    channelSelector->activateButtons();
}

Visualizer* PhaseCalculatorEditor::createNewCanvas()
{
    canvas = new PhaseCalculatorCanvas(static_cast<PhaseCalculator*>(getProcessor()));
    return canvas;
}

void PhaseCalculatorEditor::updateSettings()
{
    updateVisualizer();
}

void PhaseCalculatorEditor::saveCustomParameters(XmlElement* xml)
{
    VisualizerEditor::saveCustomParameters(xml);

    xml->setAttribute("Type", "PhaseCalculatorEditor");
    PhaseCalculator* processor = (PhaseCalculator*)(getProcessor());
    
    XmlElement* paramValues = xml->createNewChildElement("VALUES");
    paramValues->setAttribute("hilbertLength", processor->hilbertLength);
    paramValues->setAttribute("predLength", processor->predictionLength);
    paramValues->setAttribute("calcInterval", processor->calcInterval);
    paramValues->setAttribute("arOrder", processor->arOrder);
    paramValues->setAttribute("lowCut", processor->lowCut);
    paramValues->setAttribute("highCut", processor->highCut);
    paramValues->setAttribute("outputMode", processor->outputMode);
}

void PhaseCalculatorEditor::loadCustomParameters(XmlElement* xml)
{
    VisualizerEditor::loadCustomParameters(xml);

    forEachXmlChildElementWithTagName(*xml, xmlNode, "VALUES")
    {
        // some parameters have two fallbacks for backwards compatability
        hilbertLengthBox->setText(xmlNode->getStringAttribute("hilbertLength", 
            xmlNode->getStringAttribute("processLength", hilbertLengthBox->getText())), sendNotificationSync);
        predLengthEditable->setText(xmlNode->getStringAttribute("predLength",
            xmlNode->getStringAttribute("numFuture", predLengthEditable->getText())), sendNotificationSync);
        recalcIntervalEditable->setText(xmlNode->getStringAttribute("calcInterval", recalcIntervalEditable->getText()), sendNotificationSync);
        arOrderEditable->setText(xmlNode->getStringAttribute("arOrder", arOrderEditable->getText()), sendNotificationSync);
        lowCutEditable->setText(xmlNode->getStringAttribute("lowCut", lowCutEditable->getText()), sendNotificationSync);
        highCutEditable->setText(xmlNode->getStringAttribute("highCut", highCutEditable->getText()), sendNotificationSync);
        outputModeBox->setSelectedId(xmlNode->getIntAttribute("outputMode", outputModeBox->getSelectedId()), sendNotificationSync);
    }
}

void PhaseCalculatorEditor::setHighCut(float newHighCut)
{
    highCutEditable->setText(String(newHighCut), sendNotificationSync);
}

void PhaseCalculatorEditor::setVisContinuousChan(int chan)
{
    if (canvas != nullptr)
    {
        static_cast<PhaseCalculatorCanvas*>(canvas.get())->setContinuousChannel(chan);
    }
}

// static utilities

/* Attempts to parse the current text of a label as an int between min and max inclusive.
*  If successful, sets "*out" and the label text to this value and and returns true.
*  Otherwise, sets the label text to defaultValue and returns false.
*/
bool PhaseCalculatorEditor::updateIntLabel(Label* label, int min, int max, int defaultValue, int* out)
{
    const String& in = label->getText();
    int parsedInt;
    try
    {
        parsedInt = std::stoi(in.toRawUTF8());
    }
    catch (const std::logic_error&)
    {
        label->setText(String(defaultValue), dontSendNotification);
        return false;
    }

    *out = jmax(min, jmin(max, parsedInt));

    label->setText(String(*out), dontSendNotification);
    return true;
}

// Like updateIntLabel, but for floats
bool PhaseCalculatorEditor::updateFloatLabel(Label* label, float min, float max,
    float defaultValue, float* out)
{
    const String& in = label->getText();
    float parsedFloat;
    try
    {
        parsedFloat = std::stof(in.toRawUTF8());
    }
    catch (const std::logic_error&)
    {
        label->setText(String(defaultValue), dontSendNotification);
        return false;
    }

    *out = jmax(min, jmin(max, parsedFloat));

    label->setText(String(*out), dontSendNotification);
    return true;
}

// ProcessBufferSlider definitions
ProcessBufferSlider::ProcessBufferSlider(const String& componentName)
    : Slider        (componentName)
    , realMinValue  (0)
{
    setLookAndFeel(&myLookAndFeel);
    setSliderStyle(LinearBar);
    setTextBoxStyle(NoTextBox, false, 40, 20);
    setScrollWheelEnabled(false);
}

ProcessBufferSlider::~ProcessBufferSlider() {}

double ProcessBufferSlider::snapValue(double attemptedValue, DragMode dragMode)
{
    return jmax(attemptedValue, realMinValue);
}

void ProcessBufferSlider::updateFromProcessor(PhaseCalculator* parentNode)
{
    int hilbertLength = parentNode->hilbertLength;
    int predLength = parentNode->predictionLength;
    realMinValue = parentNode->arOrder;

    setRange(0, hilbertLength, 1);

    // enforce min value via snapValue and slider listener
    setValue(0); // hack to ensure the listener gets called even if only the range is changed
    setValue(hilbertLength - predLength);
}

double ProcessBufferSlider::getRealMinValue()
{
    return realMinValue;
}