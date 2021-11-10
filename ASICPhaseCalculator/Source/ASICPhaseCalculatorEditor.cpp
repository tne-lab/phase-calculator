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

#include "ASICPhaseCalculatorEditor.h"
#include "ASICPhaseCalculatorCanvas.h"
#include "HTransformers.h"
#include <climits> // INT_MAX
#include <cfloat>  // FLT_MAX
#include <cmath>   // abs

namespace ASICPhaseCalculator
{ //added by sumedh
    Editor::Editor(Node* parentNode, bool useDefaultParameterEditors)
        : VisualizerEditor(parentNode, 250, useDefaultParameterEditors)
        , extraChanManager(parentNode)
        , prevExtraChans(0)
    {
        tabText = "Event Phase Plot";
        int filterWidth = 5;

        // make the canvas now, so that restoring its parameters always works.
        canvas = new Canvas(parentNode);

        bandLabel = new Label("bandL", "Freq range:");
        bandLabel->setBounds(5, 25, 100, 20);
        bandLabel->setFont({ "Small Text", 12, Font::plain });
        bandLabel->setColour(Label::textColourId, Colours::darkgrey);
        addAndMakeVisible(bandLabel);

        bandBox = new ComboBox("bandB");
        for (int b = 0; b < NUM_BANDS; ++b)
        {
            bandBox->addItem(Hilbert::bandName[b], b + 1);
        }
        bandBox->setSelectedId(parentNode->band + 1);
        bandBox->setTooltip(freqRangeTooltip);
        bandBox->setBounds(7, 42, 110, 20);
        bandBox->addListener(this);
        addAndMakeVisible(bandBox);

        outputModeLabel = new Label("outputModeL", "Output:");
        outputModeLabel->setBounds(filterWidth, 72, 70, 20);
        outputModeLabel->setFont({ "Small Text", 12, Font::plain });
        outputModeLabel->setColour(Label::textColourId, Colours::darkgrey);
        addAndMakeVisible(outputModeLabel);

        outputModeBox = new ComboBox("outputModeB");
        outputModeBox->addItem("PHASE", PH);
        outputModeBox->addItem("MAG", MAG);
        outputModeBox->addItem("PH+MAG", PH_AND_MAG);
        outputModeBox->addItem("IMAG", IM);
        outputModeBox->setSelectedId(parentNode->outputMode);
        outputModeBox->setTooltip(outputModeTooltip);
        outputModeBox->setBounds(filterWidth + 5, 92, 76, 19);
        outputModeBox->addListener(this);
        addAndMakeVisible(outputModeBox);

        // new channels should be disabled by default
        channelSelector->paramButtonsToggledByDefault(false);
    }

    Editor::~Editor() {}

    void Editor::comboBoxChanged(ComboBox* comboBoxThatHasChanged)
    {
        Node* processor = static_cast<Node*>(getProcessor());

        if (comboBoxThatHasChanged == bandBox)
        {
            processor->setParameter(BAND, static_cast<float>(bandBox->getSelectedId() - 1));
        }
        else if (comboBoxThatHasChanged == outputModeBox)
        {
            processor->setParameter(OUTPUT_MODE, static_cast<float>(outputModeBox->getSelectedId()));
        }
    }

    void Editor::labelTextChanged(Label* labelThatHasChanged)
    {
        Node* processor = static_cast<Node*>(getProcessor());
    }

    void Editor::channelChanged(int chan, bool newState)
    {
        auto pc = static_cast<Node*>(getProcessor());
        if (chan < pc->getNumInputs())
        {
            // mark channel as activated or deactivated
            if (newState)
            {
                // check whether channel can be activated
                if (!pc->activateInputChannel(chan))
                {
                    pc->deselectChannel(chan, true);
                    return;
                }
            }
            else
            {
                pc->deactivateInputChannel(chan);
            }
            
            // update visualizer and maybe extra output channels
            if (pc->outputMode == PH_AND_MAG)
            {
                if (newState)
                {
                    extraChanManager.addExtraChan(chan);
                }
                else
                {
                    extraChanManager.removeExtraChan(chan);
                }

                // Update signal chain to add/remove output channels if necessary
                CoreServices::updateSignalChain(this);
            }
            else // (if not updating the whole signal chain)
            {
                // update the available continuous channels for visualizer
                updateVisualizer();
            }
        }
    }

    void Editor::startAcquisition()
    {
        bandBox->setEnabled(true);
        outputModeBox->setEnabled(false);
        channelSelector->inactivateButtons();
    }

    void Editor::stopAcquisition()
    {
        bandBox->setEnabled(true);
        outputModeBox->setEnabled(true);
        channelSelector->activateButtons();
    }

    Visualizer* Editor::createNewCanvas()
    {
        return canvas;
    }

    void Editor::updateSettings()
    {
        auto pc = static_cast<Node*>(getProcessor());

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

    void Editor::saveCustomParameters(XmlElement* xml)
    {
        VisualizerEditor::saveCustomParameters(xml);

        xml->setAttribute("Type", "ASICPhaseCalculatorEditor");
        Node* processor = (Node*)(getProcessor());

        XmlElement* paramValues = xml->createNewChildElement("VALUES");
        paramValues->setAttribute("lowCut", processor->lowCut);
        paramValues->setAttribute("highCut", processor->highCut);
        paramValues->setAttribute("outputMode", processor->outputMode);

        const Array<float>& validBand = Hilbert::validBand[processor->band];
        paramValues->setAttribute("rangeMin", validBand[0]);
        paramValues->setAttribute("rangeMax", validBand[1]);
    }

    void Editor::loadCustomParameters(XmlElement* xml)
    {
        VisualizerEditor::loadCustomParameters(xml);

        forEachXmlChildElementWithTagName(*xml, xmlNode, "VALUES")
        {
            // some parameters have two fallbacks for backwards compatability
            bandBox->setSelectedId(selectBandFromSavedParams(xmlNode) + 1, sendNotificationSync);
            outputModeBox->setSelectedId(xmlNode->getIntAttribute("outputMode", outputModeBox->getSelectedId()), sendNotificationSync);
        }
    }


    void Editor::refreshVisContinuousChan()
    {
        auto p = static_cast<Node*>(getProcessor());
        if (canvas != nullptr)
        {
            auto c = static_cast<Canvas*>(canvas.get());
            c->displayContinuousChan(p->visContinuousChannel);
        }
    }

    // static utilities

    int Editor::selectBandFromSavedParams(const XmlElement* xmlNode)
    {
        if (xmlNode->hasAttribute("rangeMin") && xmlNode->hasAttribute("rangeMax"))
        {
            // I don't trust JUCE's string parsing functionality (doesn't indicate failure, etc...)
            float rangeMin, rangeMax;
            if (readNumber(xmlNode->getStringAttribute("rangeMin"), rangeMin) &&
                readNumber(xmlNode->getStringAttribute("rangeMax"), rangeMax))
            {
                // try to find a matching band.
                for (int band = 0; band < NUM_BANDS; ++band)
                {
                    const Array<float>& validRange = Hilbert::validBand[band];

                    if (validRange[0] == rangeMin && validRange[1] == rangeMax)
                    {
                        return band;
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

                for (int band = 0; band < NUM_BANDS; ++band)
                {
                    const Array<float>& validRange = Hilbert::validBand[band];

                    float midBand = (validRange[0] + validRange[1]) / 2;
                    float midDist = std::abs(midBand - midCut);
                    if (validRange[0] <= lowCut && validRange[1] >= highCut && midDist < bestBandDist)
                    {
                        bestBand = band;
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

    Editor::ExtraChanManager::ExtraChanManager(const Node* processor)
        : p(processor)
    {}

    void Editor::ExtraChanManager::buttonClicked(Button* button)
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

    void Editor::ExtraChanManager::addExtraChan(int inputChan)
    {
        Array<int> activeInputs = p->getActiveInputs();
        int newInputIndex = activeInputs.indexOf(inputChan);
        jassert(newInputIndex <= recordStatus.size());
        recordStatus.insert(newInputIndex, false);
    }

    void Editor::ExtraChanManager::removeExtraChan(int inputChan)
    {
        // find # of lower-index active channels
        Array<int> activeInputs = p->getActiveInputs();
        int numActiveInputs = activeInputs.size();
        int i = 0;
        for (; i < numActiveInputs && activeInputs[i] < inputChan; ++i);
        jassert(i < recordStatus.size());
        recordStatus.remove(i);
    }

    void Editor::ExtraChanManager::resize(int numExtraChans)
    {
        recordStatus.resize(numExtraChans);
    }

    bool Editor::ExtraChanManager::getRecordStatus(int extraChan) const
    {
        return recordStatus[extraChan];
    }
}