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
#include "HTransformers.h"
#include <climits> // INT_MAX
#include <cfloat>  // FLT_MAX
#include <cmath>   // abs

PhaseCalculatorEditor::PhaseCalculatorEditor(PhaseCalculator* parentNode, bool useDefaultParameterEditors)
    : VisualizerEditor  (parentNode, 220, useDefaultParameterEditors)
    , extraChanManager  (parentNode)
    , prevExtraChans    (0)
{
    tabText = "Event Phase Plot";
    int filterWidth = 120;

    // make the canvas now, so that restoring its parameters always works.
    canvas = new PhaseCalculatorCanvas(parentNode);

    bandLabel = new Label("bandL", "Freq range:");
    bandLabel->setBounds(5, 25, 100, 20);
    bandLabel->setFont({ "Small Text", 12, Font::plain });
    bandLabel->setColour(Label::textColourId, Colours::darkgrey);
    addAndMakeVisible(bandLabel);

    bandBox = new ComboBox("bandB");
    for (int b = 0; b < NUM_BANDS; ++b)
    {
        bandBox->addItem(Hilbert::BAND_NAME.at(b), b + 1);
    }
    bandBox->setSelectedId(parentNode->band + 1);
    bandBox->setTooltip(FREQ_RANGE_TOOLTIP);
    bandBox->setBounds(7, 42, 110, 20);
    bandBox->addListener(this);
    addAndMakeVisible(bandBox);

    lowCutLabel = new Label("lowCutL", "Low:");
    lowCutLabel->setBounds(5, 70, 50, 20);
    lowCutLabel->setFont({ "Small Text", 12, Font::plain });
    lowCutLabel->setColour(Label::textColourId, Colours::darkgrey);
    addAndMakeVisible(lowCutLabel);

    lowCutEditable = new Label("lowCutE");
    lowCutEditable->setEditable(true);
    lowCutEditable->addListener(this);
    lowCutEditable->setBounds(50, 70, 35, 18);
    lowCutEditable->setText(String(parentNode->lowCut), dontSendNotification);
    lowCutEditable->setColour(Label::backgroundColourId, Colours::grey);
    lowCutEditable->setColour(Label::textColourId, Colours::white);
    addAndMakeVisible(lowCutEditable);

    lowCutUnit = new Label("lowCutU", "Hz");
    lowCutUnit->setBounds(85, 70, 25, 18);
    lowCutUnit->setFont({ "Small Text", 12, Font::plain });
    lowCutUnit->setColour(Label::textColourId, Colours::darkgrey);
    addAndMakeVisible(lowCutUnit);

    highCutLabel = new Label("highCutL", "High:");
    highCutLabel->setBounds(5, 100, 50, 20);
    highCutLabel->setFont({ "Small Text", 12, Font::plain });
    highCutLabel->setColour(Label::textColourId, Colours::darkgrey);
    addAndMakeVisible(highCutLabel);

    highCutEditable = new Label("highCutE");
    highCutEditable->setEditable(true);
    highCutEditable->addListener(this);
    highCutEditable->setBounds(50, 100, 35, 18);
    highCutEditable->setText(String(parentNode->highCut), dontSendNotification);
    highCutEditable->setColour(Label::backgroundColourId, Colours::grey);
    highCutEditable->setColour(Label::textColourId, Colours::white);
    addAndMakeVisible(highCutEditable);

    highCutUnit = new Label("highCutU", "Hz");
    highCutUnit->setBounds(85, 100, 25, 18);
    highCutUnit->setFont({ "Small Text", 12, Font::plain });
    highCutUnit->setColour(Label::textColourId, Colours::darkgrey);
    addAndMakeVisible(highCutUnit);

    recalcIntervalLabel = new Label("recalcL", "AR Refresh:");
    recalcIntervalLabel->setBounds(filterWidth, 25, 100, 20);
    recalcIntervalLabel->setFont({ "Small Text", 12, Font::plain });
    recalcIntervalLabel->setColour(Label::textColourId, Colours::darkgrey);
    addAndMakeVisible(recalcIntervalLabel);

    recalcIntervalEditable = new Label("recalcE");
    recalcIntervalEditable->setEditable(true);
    recalcIntervalEditable->addListener(this);
    recalcIntervalEditable->setBounds(filterWidth + 5, 44, 55, 18);
    recalcIntervalEditable->setColour(Label::backgroundColourId, Colours::grey);
    recalcIntervalEditable->setColour(Label::textColourId, Colours::white);
    recalcIntervalEditable->setText(String(parentNode->calcInterval), dontSendNotification);
    recalcIntervalEditable->setTooltip(RECALC_INTERVAL_TOOLTIP);
    addAndMakeVisible(recalcIntervalEditable);

    recalcIntervalUnit = new Label("recalcU", "ms");
    recalcIntervalUnit->setBounds(filterWidth + 60, 47, 25, 15);
    recalcIntervalUnit->setFont({ "Small Text", 12, Font::plain });
    recalcIntervalUnit->setColour(Label::textColourId, Colours::darkgrey);
    addAndMakeVisible(recalcIntervalUnit);

    arOrderLabel = new Label("arOrderL", "Order:");
    arOrderLabel->setBounds(filterWidth, 65, 60, 20);
    arOrderLabel->setFont({ "Small Text", 12, Font::plain });
    arOrderLabel->setColour(Label::textColourId, Colours::darkgrey);
    addAndMakeVisible(arOrderLabel);

    arOrderEditable = new Label("arOrderE");
    arOrderEditable->setEditable(true);
    arOrderEditable->addListener(this);
    arOrderEditable->setBounds(filterWidth + 55, 66, 25, 18);
    arOrderEditable->setColour(Label::backgroundColourId, Colours::grey);
    arOrderEditable->setColour(Label::textColourId, Colours::white);
    arOrderEditable->setText(String(parentNode->arOrder), sendNotificationAsync);
    arOrderEditable->setTooltip(AR_ORDER_TOOLTIP);
    addAndMakeVisible(arOrderEditable);

    outputModeLabel = new Label("outputModeL", "Output:");
    outputModeLabel->setBounds(filterWidth, 87, 70, 20);
    outputModeLabel->setFont({ "Small Text", 12, Font::plain });
    outputModeLabel->setColour(Label::textColourId, Colours::darkgrey);
    addAndMakeVisible(outputModeLabel);

    outputModeBox = new ComboBox("outputModeB");
    outputModeBox->addItem("PHASE", PH);
    outputModeBox->addItem("MAG", MAG);
    outputModeBox->addItem("PH+MAG", PH_AND_MAG);
    outputModeBox->addItem("IMAG", IM);
    outputModeBox->setSelectedId(parentNode->outputMode);
    outputModeBox->setTooltip(OUTPUT_MODE_TOOLTIP);
    outputModeBox->setBounds(filterWidth + 5, 105, 76, 19);
    outputModeBox->addListener(this);
    addAndMakeVisible(outputModeBox);

    // new channels should be disabled by default
    channelSelector->paramButtonsToggledByDefault(false);
}

PhaseCalculatorEditor::~PhaseCalculatorEditor() {}

void PhaseCalculatorEditor::comboBoxChanged(ComboBox* comboBoxThatHasChanged)
{
    PhaseCalculator* processor = static_cast<PhaseCalculator*>(getProcessor());

    if (comboBoxThatHasChanged == bandBox)
    {
        processor->setParameter(BAND, static_cast<float>(bandBox->getSelectedId() - 1));
    }
    else if (comboBoxThatHasChanged == outputModeBox)
    {
        processor->setParameter(OUTPUT_MODE, static_cast<float>(outputModeBox->getSelectedId()));
    }
}

void PhaseCalculatorEditor::labelTextChanged(Label* labelThatHasChanged)
{
    PhaseCalculator* processor = static_cast<PhaseCalculator*>(getProcessor());

    if (labelThatHasChanged == recalcIntervalEditable)
    {
        int intInput;
        bool valid = updateControl(labelThatHasChanged, 0, INT_MAX, processor->calcInterval, intInput);

        if (valid)
        {
            processor->setParameter(RECALC_INTERVAL, static_cast<float>(intInput));
        }
    }
    else if (labelThatHasChanged == arOrderEditable)
    {
        int intInput;
        bool valid = updateControl(labelThatHasChanged, 1, INT_MAX, processor->arOrder, intInput);

        if (valid)
        {
            processor->setParameter(AR_ORDER, static_cast<float>(intInput));
        }
    }
    else if (labelThatHasChanged == lowCutEditable)
    {
        float floatInput;
        bool valid = updateControl(labelThatHasChanged, 0.0f, FLT_MAX, processor->lowCut, floatInput);

        if (valid)
        {
            processor->setParameter(LOWCUT, floatInput);
        }
    }
    else if (labelThatHasChanged == highCutEditable)
    {
        float floatInput;
        bool valid = updateControl(labelThatHasChanged, 0.0f, FLT_MAX, processor->highCut, floatInput);
        
        if (valid)
        {
            processor->setParameter(HIGHCUT, floatInput);
        }
    }
}

void PhaseCalculatorEditor::channelChanged(int chan, bool newState)
{
    auto pc = static_cast<PhaseCalculator*>(getProcessor());
    if (chan < pc->getNumInputs())
    {
        Array<int> activeInputs = pc->getActiveInputs();
        if (newState)
        {
            // check whether sample rate is compatible (and if not, disable channel)
            if (!pc->validateSampleRate(chan)) { return; }

            // ensure space allocated for per-active-channel arrays
            if (activeInputs.size() > pc->numActiveChansAllocated)
            {
                pc->addActiveChannel();
            }
        }

        if (pc->outputMode == PH_AND_MAG)
        {
            if (newState)
            {
                extraChanManager.addExtraChan(chan, activeInputs);
            }
            else
            {
                extraChanManager.removeExtraChan(chan, activeInputs);
            }

            // Update signal chain to add/remove output channels if necessary
            CoreServices::updateSignalChain(this);
        }
        else
        {
            updateVisualizer(); // update the available continuous channels for visualizer
        }
    }
}

void PhaseCalculatorEditor::startAcquisition()
{
    bandBox->setEnabled(false);
    lowCutEditable->setEnabled(false);
    highCutEditable->setEnabled(false);
    arOrderEditable->setEnabled(false);
    outputModeBox->setEnabled(false);
    channelSelector->inactivateButtons();
}

void PhaseCalculatorEditor::stopAcquisition()
{
    bandBox->setEnabled(true);
    lowCutEditable->setEnabled(true);
    highCutEditable->setEnabled(true);
    arOrderEditable->setEnabled(true);
    outputModeBox->setEnabled(true);
    channelSelector->activateButtons();
}

Visualizer* PhaseCalculatorEditor::createNewCanvas()
{
    return canvas;
}

void PhaseCalculatorEditor::updateSettings()
{
    auto pc = static_cast<PhaseCalculator*>(getProcessor());

    // only care about any of this stuff if we have extra channels
    // (and preserve when deselecting/reselecting PH_AND_MAG)
    if (pc->outputMode != PH_AND_MAG || channelSelector == nullptr) { return; }
    
    int numChans = pc->getNumOutputs();
    int numInputs = pc->getNumInputs();
    int extraChans = numChans - numInputs;

    int prevNumChans = channelSelector->getNumChannels();
    int prevNumInputs = prevNumChans - prevExtraChans;
    prevExtraChans = extraChans; // update for next time

    extraChanManager.resize(extraChans);
    channelSelector->setNumChannels(numChans);

    // super hacky, access record buttons to add or remove listeners
    Component* rbmComponent = channelSelector->getChildComponent(9);
    auto recordButtonManager = dynamic_cast<ButtonGroupManager*>(rbmComponent);
    if (recordButtonManager == nullptr)
    {
        jassertfalse;
        return;
    }

    // remove listeners on channels that are no longer "extra channels"
    // and set their record status to false since they're actually new channels
    for (int chan = prevNumInputs; chan < jmin(prevNumChans, numInputs); ++chan)
    {
        juce::Button* recordButton = recordButtonManager->getButtonAt(chan);
        recordButton->removeListener(&extraChanManager);
        // make sure listener really gets called
        recordButton->setToggleState(true, dontSendNotification);
        channelSelector->setRecordStatus(chan, false);
    }

    // add listeners for current "extra channels" and restore record statuses
    // (it's OK if addListener gets called more than once for a button)
    for (int eChan = 0; eChan < extraChans; ++eChan)
    {
        int chan = numInputs + eChan;
        juce::Button* recordButton = recordButtonManager->getButtonAt(chan);
        recordButton->removeListener(&extraChanManager);
        // make sure listener really gets called
        bool recordStatus = extraChanManager.getRecordStatus(eChan);
        recordButton->setToggleState(!recordStatus, dontSendNotification);
        channelSelector->setRecordStatus(chan, recordStatus);
        recordButton->addListener(&extraChanManager);
    }
}

void PhaseCalculatorEditor::saveCustomParameters(XmlElement* xml)
{
    VisualizerEditor::saveCustomParameters(xml);

    xml->setAttribute("Type", "PhaseCalculatorEditor");
    PhaseCalculator* processor = (PhaseCalculator*)(getProcessor());

    XmlElement* paramValues = xml->createNewChildElement("VALUES");
    paramValues->setAttribute("calcInterval", processor->calcInterval);
    paramValues->setAttribute("arOrder", processor->arOrder);
    paramValues->setAttribute("lowCut", processor->lowCut);
    paramValues->setAttribute("highCut", processor->highCut);
    paramValues->setAttribute("outputMode", processor->outputMode);

    const Array<float>& validBand = Hilbert::VALID_BAND.at(processor->band);
    paramValues->setAttribute("rangeMin", validBand[0]);
    paramValues->setAttribute("rangeMax", validBand[1]);
}

void PhaseCalculatorEditor::loadCustomParameters(XmlElement* xml)
{
    VisualizerEditor::loadCustomParameters(xml);

    forEachXmlChildElementWithTagName(*xml, xmlNode, "VALUES")
    {
        // some parameters have two fallbacks for backwards compatability
        recalcIntervalEditable->setText(xmlNode->getStringAttribute("calcInterval", recalcIntervalEditable->getText()), sendNotificationSync);
        arOrderEditable->setText(xmlNode->getStringAttribute("arOrder", arOrderEditable->getText()), sendNotificationSync);
        bandBox->setSelectedId(selectBandFromSavedParams(xmlNode) + 1, sendNotificationSync);
        lowCutEditable->setText(xmlNode->getStringAttribute("lowCut", lowCutEditable->getText()), sendNotificationSync);
        highCutEditable->setText(xmlNode->getStringAttribute("highCut", highCutEditable->getText()), sendNotificationSync);
        outputModeBox->setSelectedId(xmlNode->getIntAttribute("outputMode", outputModeBox->getSelectedId()), sendNotificationSync);
    }
}

void PhaseCalculatorEditor::refreshLowCut()
{
    auto p = static_cast<PhaseCalculator*>(getProcessor());
    lowCutEditable->setText(String(p->lowCut), dontSendNotification);
}

void PhaseCalculatorEditor::refreshHighCut()
{
    auto p = static_cast<PhaseCalculator*>(getProcessor());
    highCutEditable->setText(String(p->highCut), dontSendNotification);
}

void PhaseCalculatorEditor::refreshVisContinuousChan()
{
    auto p = static_cast<PhaseCalculator*>(getProcessor());
    if (canvas != nullptr)
    {
        auto c = static_cast<PhaseCalculatorCanvas*>(canvas.get());
        c->displayContinuousChan(p->visContinuousChannel);
    }
}

// static utilities

int PhaseCalculatorEditor::selectBandFromSavedParams(const XmlElement* xmlNode)
{
    if (xmlNode->hasAttribute("rangeMin") && xmlNode->hasAttribute("rangeMax"))
    {
        // I don't trust JUCE's string parsing functionality (doesn't indicate failure, etc...)
        float rangeMin, rangeMax;
        if (readNumber(xmlNode->getStringAttribute("rangeMin"), rangeMin) &&
            readNumber(xmlNode->getStringAttribute("rangeMax"), rangeMax))
        {
            // try to find a matching band.
            for (auto keyVal : Hilbert::VALID_BAND)
            {
                int b = keyVal.first;
                const Array<float>& validRange = keyVal.second;

                if (validRange[0] == rangeMin && validRange[1] == rangeMax)
                {
                    return b;
                }
            }
        }
    }

    // no match, or rangeMin and rangeMax were invalid
    // try to find a good fit for lowCut and highCut
    if (xmlNode->hasAttribute("lowCut") && xmlNode->hasAttribute("highCut"))
    {
        float lowCut, highCut;
        if (readNumber(xmlNode->getStringAttribute("lowCut"), lowCut) && lowCut >= 0 &&
            readNumber(xmlNode->getStringAttribute("highCut"), highCut) && highCut > lowCut)
        {
            float midCut = (lowCut + highCut) / 2;
            int bestBand = 0;
            float bestBandDist = FLT_MAX;
            
            for (auto keyVal : Hilbert::VALID_BAND)
            {
                int b = keyVal.first;
                const Array<float>& validRange = keyVal.second;

                float midBand = (validRange[0] + validRange[1]) / 2;
                float midDist = std::abs(midBand - midCut);
                if (validRange[0] <= lowCut && validRange[1] >= highCut && midDist < bestBandDist)
                {
                    bestBand = b;
                    bestBandDist = midDist;
                }
            }
            return bestBand;
        }
    }

    // just fallback to first band
    return 0;
}

// -------- ExtraChanManager ---------

PhaseCalculatorEditor::ExtraChanManager::ExtraChanManager(const PhaseCalculator* processor)
    : p(processor)
{}

void PhaseCalculatorEditor::ExtraChanManager::buttonClicked(Button* button)
{
    int numInputs = p->getNumInputs();
    int chanInd = button->getParentComponent()->getIndexOfChildComponent(button);
    int extraChanInd = chanInd - numInputs;
    if (extraChanInd < 0 || extraChanInd >= recordStatus.size())
    {
        jassertfalse;
        return;
    }
    recordStatus.set(extraChanInd, button->getToggleState());
}

void PhaseCalculatorEditor::ExtraChanManager::addExtraChan(int inputChan, const Array<int>& activeInputs)
{
    int newInputIndex = activeInputs.indexOf(inputChan);
    jassert(newInputIndex <= recordStatus.size());
    recordStatus.insert(newInputIndex, false);
}

void PhaseCalculatorEditor::ExtraChanManager::removeExtraChan(int inputChan, const Array<int>& activeInputs)
{
    // find # of lower-index active inputs
    int i = 0;
    int numActiveInputs = activeInputs.size();
    for (; i < numActiveInputs && activeInputs[i] < inputChan; ++i);
    jassert(i < recordStatus.size());
    recordStatus.remove(i);
}

void PhaseCalculatorEditor::ExtraChanManager::resize(int numExtraChans)
{
    recordStatus.resize(numExtraChans);
}

bool PhaseCalculatorEditor::ExtraChanManager::getRecordStatus(int extraChan) const
{
    return recordStatus[extraChan];
}