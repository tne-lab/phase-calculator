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

#include "AnalyticSignalEditor.h"
#include "AnalyticSignalCanvas.h"
#include "HTransformers.h"
#include <climits> // INT_MAX
#include <cfloat>  // FLT_MAX
#include <cmath>   // abs

namespace StatePhaseEst
{
    Editor::Editor(Node* parentNode, bool useDefaultParameterEditors)
        : VisualizerEditor(parentNode, 220, useDefaultParameterEditors = false)
        , extraChanManager(parentNode)
        , prevExtraChans(0)
    {

        tabText = "Event Phase Plot";

        // make the canvas now, so that restoring its parameters always works.
        canvas = new Canvas(parentNode);

        switch (parentNode->curPhaseAlg)
        {
        case UNKNOWN_PHASE_ALG:
            createSelectPhaseAlgComponents();
            break;

        case HILBERT_TRANSFORMER:
            createHilbertComponents();
            break;

        case STATE_SPACE:
            createSSPEComponents();
            break;

        }


        // new channels should be disabled by default
        channelSelector->paramButtonsToggledByDefault(false);

    }

    Editor::~Editor() {}

    void Editor::createSelectPhaseAlgComponents()
    {
        setDesiredWidth(220);
        int filterWidth = 120;
        auto p = static_cast<Node*>(getProcessor());

        phaseSelectLabel = new Label("phaseSelectL", "Select which phase algorithm to use:");
        phaseSelectLabel->setBounds(5, 25, 200, 20);
        phaseSelectLabel->setFont({ "Small Text", 12, Font::plain });
        phaseSelectLabel->setColour(Label::textColourId, Colours::darkgrey);
        addAndMakeVisible(phaseSelectLabel);

        phaseSelectBox = new ComboBox("phaseSelectB");
        phaseSelectBox->addItem("HILBERT_TRANSFORMER", 1);
        phaseSelectBox->addItem("STATE_SPACE", 2);
        phaseSelectBox->setSelectedId(p->curPhaseAlg);
        phaseSelectBox->setTooltip(freqRangeTooltip);
        phaseSelectBox->setBounds(7, 42, 200, 20);
        phaseSelectBox->addListener(this);
        addAndMakeVisible(phaseSelectBox);
    }


    void Editor::createHilbertComponents()
    {

        setDesiredWidth(220);
        int filterWidth = 120;
        auto p = static_cast<Node*>(getProcessor());

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
        bandBox->setSelectedId(p->band + 1);
        bandBox->setTooltip(freqRangeTooltip);
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
        lowCutEditable->setText(String(p->lowCut), dontSendNotification);
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
        highCutEditable->setText(String(p->highCut), dontSendNotification);
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
        recalcIntervalEditable->setText(String(p->calcInterval), dontSendNotification);
        recalcIntervalEditable->setTooltip(recalcIntervalTooltip);
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
        arOrderEditable->setText(String(p->arOrder), sendNotificationAsync);
        arOrderEditable->setTooltip(arOrderTooltip);
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
        outputModeBox->setSelectedId(p->outputMode);
        outputModeBox->setTooltip(outputModeTooltip);
        outputModeBox->setBounds(filterWidth + 5, 105, 76, 19);
        outputModeBox->addListener(this);
        addAndMakeVisible(outputModeBox);

    }

    void Editor::createSSPEComponents()
    {
        setDesiredWidth(220);
        
        int filterWidth = 96;
        auto p = static_cast<Node*>(getProcessor());

        numFreqsLabel = new Label("nFreqL", "Num Freqs:");
        numFreqsLabel->setBounds(5, 25, 75, 20);
        numFreqsLabel->setFont({ "Small Text", 12, Font::plain });
        numFreqsLabel->setColour(Label::textColourId, Colours::darkgrey);
        addAndMakeVisible(numFreqsLabel);

        numFreqsBox = new ComboBox("bandB");
        for (int b = 1; b <= 3; ++b)
        {
            numFreqsBox->addItem(String(b), b);
        }
        numFreqsBox->setSelectedId(p->getFreqs().size());
        //numFreqsBox->setTooltip(freqRangeTooltip);           // Change me!!!!!
        numFreqsBox->setBounds(7, 42, 68, 20);
        numFreqsBox->addListener(this);
        addAndMakeVisible(numFreqsBox);

        createFreqArrayLabelEdit(freqOneLabel, freqOneEditable, covOneLabel, covOneEditable, 0);
        createFreqArrayLabelEdit(freqTwoLabel, freqTwoEditable, covTwoLabel, covTwoEditable, 1);
        createFreqArrayLabelEdit(freqThreeLabel, freqThreeEditable, covThreeLabel, covThreeEditable, 2);


        foiLabel = new Label("foiL", "FOI");
        foiLabel->setBounds(filterWidth - 21, 50, 30, 20);
        foiLabel->setFont({ "Small Text", 12, Font::plain });
        foiLabel->setColour(Label::textColourId, Colours::darkgrey);
        addAndMakeVisible(foiLabel);

        int num = 0;
        freqOneButton = new ToggleButton();
        freqOneButton->setBounds(filterWidth - 20, 59 + 21 * num, 25, 30);
        freqOneButton->addListener(this);
        freqOneButton->setToggleState(false, dontSendNotification);
        freqOneButton->setColour(ToggleButton::textColourId, Colours::white);
        //averageButton->setTooltip("Average event related potentials");
        freqOneButton->setRadioGroupId(1, dontSendNotification);

        num++;
        freqTwoButton = new ToggleButton();
        freqTwoButton->setBounds(filterWidth - 20, 59 + 21 * num, 25, 30);
        freqTwoButton->addListener(this);
        freqTwoButton->setToggleState(false, dontSendNotification);
        freqTwoButton->setColour(ToggleButton::textColourId, Colours::white);
        //averageButton->setTooltip("Average event related potentials");
        freqTwoButton->setRadioGroupId(1, dontSendNotification);

        num++;
        freqThreeButton = new ToggleButton();
        freqThreeButton->setBounds(filterWidth - 20, 59 + 21 * num, 25, 30);
        freqThreeButton->addListener(this);
        freqThreeButton->setToggleState(false, dontSendNotification);
        freqThreeButton->setColour(ToggleButton::textColourId, Colours::white);
        //averageButton->setTooltip("Average event related potentials");
        freqThreeButton->setRadioGroupId(1, dontSendNotification);

        freqOneEditable->setText(String(p->freqs[0]), dontSendNotification);
        addAndMakeVisible(freqOneEditable);
        addAndMakeVisible(freqOneLabel);
        addAndMakeVisible(freqOneButton);

        covOneEditable->setText(String(p->qEst[0]), dontSendNotification);
        addAndMakeVisible(covOneEditable);
        addAndMakeVisible(covOneLabel);


        obsErrorLabel = new Label("obsErrorL", "Obs Err:");
        obsErrorLabel->setBounds(filterWidth, 25, 65, 20);
        obsErrorLabel->setFont({ "Small Text", 12, Font::plain });
        obsErrorLabel->setColour(Label::textColourId, Colours::darkgrey);
        addAndMakeVisible(obsErrorLabel);

        obsErrorEditable = new Label("obsErrorE");
        obsErrorEditable->setEditable(true);
        obsErrorEditable->addListener(this);
        obsErrorEditable->setText(String(p->obsErrorEst), dontSendNotification);
        obsErrorEditable->setBounds(filterWidth + 62, 25, 45, 18);
        obsErrorEditable->setColour(Label::backgroundColourId, Colours::grey);
        obsErrorEditable->setColour(Label::textColourId, Colours::white);
        addAndMakeVisible(obsErrorEditable);
        /*
        covarianceLabel = new Label("covarianceL", "Covariance:");
        covarianceLabel->setBounds(filterWidth, 25 + 21, 65, 20);
        covarianceLabel->setFont({ "Small Text", 12, Font::plain });
        covarianceLabel->setColour(Label::textColourId, Colours::darkgrey);
        addAndMakeVisible(covarianceLabel);

        covarianceEditable = new Label("covarianceE");
        covarianceEditable->setEditable(true);
        covarianceEditable->addListener(this);
        covarianceEditable->setText(String(p->obsErrorEst), dontSendNotification);
        covarianceEditable->setBounds(filterWidth + 62, 25 + 21, 45, 18);
        covarianceEditable->setColour(Label::backgroundColourId, Colours::grey);
        covarianceEditable->setColour(Label::textColourId, Colours::white);
        addAndMakeVisible(covarianceEditable);
        */
        // Set defualt vals
        freqOneButton->setToggleState(false, dontSendNotification);
        freqTwoButton->setToggleState(false, dontSendNotification);
        freqThreeButton->setToggleState(false, dontSendNotification);

        switch (p->freqs.size())
        {
        case 3:
            addAndMakeVisible(freqThreeLabel);
            addAndMakeVisible(freqThreeEditable);
            addAndMakeVisible(freqThreeButton);
            freqThreeEditable->setText(String(p->freqs[2]), dontSendNotification);

            addAndMakeVisible(covThreeLabel);
            addAndMakeVisible(covThreeEditable);
            covThreeEditable->setText(String(p->qEst[2]), dontSendNotification);

            // drop through
        case 2:
            addAndMakeVisible(freqTwoLabel);
            addAndMakeVisible(freqTwoEditable);
            addAndMakeVisible(freqTwoButton);
            freqTwoEditable->setText(String(p->freqs[1]), dontSendNotification);

            addAndMakeVisible(covTwoLabel);
            addAndMakeVisible(covTwoEditable);
            covTwoEditable->setText(String(p->qEst[1]), dontSendNotification);
            break;
        }

        switch (p->foi)
        {
        case 0:
            freqOneButton->setToggleState(true, dontSendNotification);
            break;
        case 1:
            freqTwoButton->setToggleState(true, dontSendNotification);
            break;
        case 2:
            freqThreeButton->setToggleState(true, dontSendNotification);
            break;
        }

    }

    void Editor::createFreqArrayLabelEdit(ScopedPointer<Label>& freqLabel, ScopedPointer<Label>& freqEditable, ScopedPointer<Label>& covLabel, ScopedPointer<Label>& covEditable, int num)
    {
        // create Freq label and editable
        freqLabel = new Label("nFreq" + String(num) + "L", "Freq:");
        freqLabel->setBounds(5, 65 + 21 * num, 50, 20);
        freqLabel->setFont({ "Small Text", 12, Font::plain });
        freqLabel->setColour(Label::textColourId, Colours::darkgrey);

        freqEditable = new Label("nFreq" + String(num) + "E");
        freqEditable->setEditable(true);
        freqEditable->addListener(this);
        freqEditable->setBounds(45, 65 + 21 * num, 33, 18);
        freqEditable->setColour(Label::backgroundColourId, Colours::grey);
        freqEditable->setColour(Label::textColourId, Colours::white);
        // freqEditable]->setTooltip(arOrderTooltip);                           // Change me!

        // create covairance label and editable
        covLabel = new Label("cov" + String(num) + "L", "Q:");
        covLabel->setBounds(115, 65 + 21 * num, 50, 20);
        covLabel->setFont({ "Small Text", 12, Font::plain });
        covLabel->setColour(Label::textColourId, Colours::darkgrey);

        covEditable = new Label("cov" + String(num) + "E");
        covEditable->setEditable(true);
        covEditable->addListener(this);
        covEditable->setBounds(135, 65 + 21 * num, 33, 18);
        covEditable->setColour(Label::backgroundColourId, Colours::grey);
        covEditable->setColour(Label::textColourId, Colours::white);

    }

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

        else if (comboBoxThatHasChanged == phaseSelectBox)
        {
            int phaseAlg = static_cast<int>(phaseSelectBox->getSelectedId());
            processor->setParameter(PHASE_ALG, phaseAlg);
            phaseSelectLabel->setVisible(false);
            phaseSelectBox->setVisible(false);
            switch (phaseAlg)
            {
            case HILBERT_TRANSFORMER:
                createHilbertComponents();
                processor->updateAllChannels();
                break;
            case STATE_SPACE:
                createSSPEComponents();
                processor->updateAllChannels();
                break;
            }
        }

        else if (comboBoxThatHasChanged == numFreqsBox)
        {
            int nFreqs = static_cast<int>(numFreqsBox->getSelectedId());
            freqThreeEditable->setVisible(false);
            freqThreeLabel->setVisible(false);
            freqThreeButton->setVisible(false);
            freqTwoEditable->setVisible(false);
            freqTwoLabel->setVisible(false);
            freqTwoButton->setVisible(false);

            covThreeEditable->setVisible(false);
            covThreeLabel->setVisible(false);
            covTwoEditable->setVisible(false);
            covTwoLabel->setVisible(false);

            processor->setParameter(N_FREQS, nFreqs);

            int foi = processor->getFoi();

            if (foi > nFreqs-1)
            {
                foi = nFreqs - 1;
                processor->setParameter(FOI, foi);
                switch (foi)
                {
                case 2:
                    freqThreeButton->setToggleState(true, sendNotificationSync);
                    break;
                case 1:
                    freqTwoButton->setToggleState(true, sendNotificationSync);
                    break;
                case 0:
                    freqOneButton->setToggleState(true, sendNotificationSync);
                    break;
                }
            }

            switch (nFreqs)
            {
            case 3:
                // freqs
                addAndMakeVisible(freqThreeLabel);
                addAndMakeVisible(freqThreeEditable);
                addAndMakeVisible(freqThreeButton);
                processor->setParameter(FREQ_THREE, freqThreeEditable->getText().getFloatValue());
                // covariance
                addAndMakeVisible(covThreeLabel);
                addAndMakeVisible(covThreeEditable);
                processor->setParameter(Q_EST_THREE, covThreeEditable->getText().getFloatValue());
                // drop through
            case 2:
                // freqs
                addAndMakeVisible(freqTwoLabel);
                addAndMakeVisible(freqTwoEditable);
                addAndMakeVisible(freqTwoButton);
                processor->setParameter(FREQ_TWO, freqTwoEditable->getText().getFloatValue());
                // covariance
                addAndMakeVisible(covTwoLabel);
                addAndMakeVisible(covTwoEditable);
                processor->setParameter(Q_EST_TWO, covTwoEditable->getText().getFloatValue());
                break;
            case 1:
                processor->updateActiveChannels(); // Just in case FOI gets updated. 
            }
        }
    }

    void Editor::labelTextChanged(Label* labelThatHasChanged)
    {
        Node* processor = static_cast<Node*>(getProcessor());

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

        else if (labelThatHasChanged == freqOneEditable)
        {
            float floatInput;
            bool valid = updateControl(labelThatHasChanged, 0.0f, FLT_MAX, processor->freqs[0], floatInput);

            if (valid)
            {
                processor->setParameter(FREQ_ONE, floatInput);
            }
        }
        else if (labelThatHasChanged == freqTwoEditable)
        {
            float floatInput;
            bool valid = updateControl(labelThatHasChanged, 0.0f, FLT_MAX, processor->freqs[1], floatInput);

            if (valid)
            {
                processor->setParameter(FREQ_TWO, floatInput);
            }
        }
        else if (labelThatHasChanged == freqThreeEditable)
        {
            float floatInput;
            bool valid = updateControl(labelThatHasChanged, 0.0f, FLT_MAX, processor->freqs[2], floatInput);

            if (valid)
            {
                processor->setParameter(FREQ_THREE, floatInput);
            }
        }
        /*
        else if (labelThatHasChanged == ampEstEditable)
        {
            float floatInput;
            bool valid = updateControl(labelThatHasChanged, 0.0f, FLT_MAX, processor->ampEst, floatInput);

            if (valid)
            {
                processor->setParameter(AMP_EST, floatInput);
            }
        }
        else if (labelThatHasChanged == qEstEditable)
        {
            float floatInput;
            bool valid = updateControl(labelThatHasChanged, 0.0f, FLT_MAX, processor->qEst, floatInput);

            if (valid)
            {
                processor->setParameter(Q_EST, floatInput);
            }
        }
        */
        else if (labelThatHasChanged == obsErrorEditable)
        {
            float floatInput;
            bool valid = updateControl(labelThatHasChanged, 0.0f, FLT_MAX, processor->obsErrorEst, floatInput);

            if (valid)
            {
                processor->setParameter(OBS_ERR_EST, floatInput);
            }
        }

        else if (labelThatHasChanged == covOneEditable)
        {
            float floatInput;
            bool valid = updateControl(labelThatHasChanged, 0.0f, FLT_MAX, processor->qEst[0], floatInput);

            if (valid)
            {
                processor->setParameter(Q_EST_ONE, floatInput);
            }
        }
        else if (labelThatHasChanged == covTwoEditable)
        {
        float floatInput;
        bool valid = updateControl(labelThatHasChanged, 0.0f, FLT_MAX, processor->qEst[1], floatInput);

        if (valid)
        {
            processor->setParameter(Q_EST_TWO, floatInput);
        }
        }
        else if (labelThatHasChanged == covThreeEditable)
        {
        float floatInput;
        bool valid = updateControl(labelThatHasChanged, 0.0f, FLT_MAX, processor->qEst[2], floatInput);

        if (valid)
        {
            processor->setParameter(Q_EST_THREE, floatInput);
        }
        }
    }

    void Editor::buttonEvent(Button* button)
    {
        Node* processor = static_cast<Node*>(getProcessor());

        if (button == freqOneButton)
        {
            processor->setParameter(FOI, 0);
            processor->updateActiveChannels();
        }

        else if (button == freqTwoButton)
        {
            processor->setParameter(FOI, 1);
            processor->updateActiveChannels();
        }
        else if (button == freqThreeButton)
        {
            processor->setParameter(FOI, 2);
            processor->updateActiveChannels();
        }
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
        auto pc = static_cast<Node*>(getProcessor());
        if (pc->curPhaseAlg == HILBERT_TRANSFORMER)
        {
            bandBox->setEnabled(false);
            lowCutEditable->setEnabled(false);
            highCutEditable->setEnabled(false);
            arOrderEditable->setEnabled(false);
            outputModeBox->setEnabled(false);
        }

        channelSelector->inactivateButtons();
    }

    void Editor::stopAcquisition()
    {
        auto pc = static_cast<Node*>(getProcessor());
        if (pc->curPhaseAlg == HILBERT_TRANSFORMER)
        {
            bandBox->setEnabled(true);
            lowCutEditable->setEnabled(true);
            highCutEditable->setEnabled(true);
            arOrderEditable->setEnabled(true);
            outputModeBox->setEnabled(true);
        }
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

        xml->setAttribute("Type", "StatePhaseEstEditor");
        Node* processor = (Node*)(getProcessor());


        XmlElement* paramValues = xml->createNewChildElement("VALUES");
        const Array<float>& validBand = Hilbert::validBand[processor->band];
        switch (processor->curPhaseAlg)
        {
        case HILBERT_TRANSFORMER:
            paramValues->setAttribute("phaseAlg", HILBERT_TRANSFORMER);
            paramValues->setAttribute("calcInterval", processor->calcInterval);
            paramValues->setAttribute("arOrder", processor->arOrder);
            paramValues->setAttribute("lowCut", processor->lowCut);
            paramValues->setAttribute("highCut", processor->highCut);
            paramValues->setAttribute("outputMode", processor->outputMode);

            
            paramValues->setAttribute("rangeMin", validBand[0]);
            paramValues->setAttribute("rangeMax", validBand[1]);
            break;
        case STATE_SPACE:
            paramValues->setAttribute("phaseAlg", STATE_SPACE);
            paramValues->setAttribute("obsErrorEst", processor->obsErrorEst);
            paramValues->setAttribute("foi", processor->getFoi());
            paramValues->setAttribute("nFreqs", processor->freqs.size());

            for (int i = 0; i < processor->freqs.size(); i++)
            {
                paramValues->setAttribute("freq_"+String(i+1), processor->freqs[i]);
                paramValues->setAttribute("qEst_" + String(i + 1), processor->qEst[i]);
               
            }
            break;
        }

        

    }

    void Editor::loadCustomParameters(XmlElement* xml)
    {
        VisualizerEditor::loadCustomParameters(xml);

        auto p = static_cast<Node*>(getProcessor());
        int cpa = 0;
        forEachXmlChildElementWithTagName(*xml, xmlNode, "VALUES")
        {
            // some parameters have two fallbacks for backwards compatability
            (xmlNode->getStringAttribute("phaseAlg", "0"), sendNotificationSync);

            cpa = xmlNode->getStringAttribute("phaseAlg", "0").getIntValue();
            p->setParameter(PHASE_ALG, cpa);

            switch (cpa)
            {
            case UNKNOWN_PHASE_ALG:
                createSelectPhaseAlgComponents();
                break;
            case HILBERT_TRANSFORMER:
                phaseSelectLabel->setVisible(false);
                phaseSelectBox->setVisible(false);
                createHilbertComponents();
               // p->setParameter(PHASE_ALG, cpa);

                recalcIntervalEditable->setText(xmlNode->getStringAttribute("calcInterval", recalcIntervalEditable->getText()), sendNotificationSync);
                arOrderEditable->setText(xmlNode->getStringAttribute("arOrder", arOrderEditable->getText()), sendNotificationSync);
                bandBox->setSelectedId(selectBandFromSavedParams(xmlNode) + 1, sendNotificationSync);
                lowCutEditable->setText(xmlNode->getStringAttribute("lowCut", lowCutEditable->getText()), sendNotificationSync);
                highCutEditable->setText(xmlNode->getStringAttribute("highCut", highCutEditable->getText()), sendNotificationSync);
                outputModeBox->setSelectedId(xmlNode->getIntAttribute("outputMode", outputModeBox->getSelectedId()), sendNotificationSync);
                break;
            case STATE_SPACE:
                phaseSelectLabel->setVisible(false);
                phaseSelectBox->setVisible(false);
                createSSPEComponents();

                obsErrorEditable->setText(xmlNode->getStringAttribute("obsErrorEst", obsErrorEditable->getText()), sendNotificationSync);
                int nFreqs = xmlNode->getStringAttribute("nFreqs", "1").getIntValue();

                numFreqsBox->setSelectedId(nFreqs, sendNotificationSync);
                switch (nFreqs)
                {
                case 3:
                    addAndMakeVisible(freqThreeLabel);
                    addAndMakeVisible(freqThreeEditable);
                    addAndMakeVisible(freqThreeButton);
                    freqThreeEditable->setText(xmlNode->getStringAttribute("freq_3", freqThreeEditable->getText()), sendNotificationSync);

                    addAndMakeVisible(covThreeEditable);
                    addAndMakeVisible(covThreeLabel);
                    covThreeEditable->setText(xmlNode->getStringAttribute("qEst_3", covThreeEditable->getText()), sendNotificationSync);
                case 2:
                    addAndMakeVisible(freqTwoLabel);
                    addAndMakeVisible(freqTwoEditable);
                    addAndMakeVisible(freqTwoButton);
                    freqTwoEditable->setText(xmlNode->getStringAttribute("freq_2", freqTwoEditable->getText()), sendNotificationSync);

                    addAndMakeVisible(covTwoEditable);
                    addAndMakeVisible(covTwoLabel);
                    covTwoEditable->setText(xmlNode->getStringAttribute("qEst_2", covTwoEditable->getText()), sendNotificationSync);
                case 1:
                    freqOneEditable->setText(xmlNode->getStringAttribute("freq_1", freqOneEditable->getText()), sendNotificationSync);
                    covOneEditable->setText(xmlNode->getStringAttribute("qEst_1", covOneEditable->getText()), sendNotificationSync);
                    break;
                }

                switch (xmlNode->getStringAttribute("foi", "0").getIntValue())
                {
                case 2:
                    freqThreeButton->setToggleState(true, sendNotificationSync);
                    break;
                case 1:
                    freqTwoButton->setToggleState(true, sendNotificationSync);
                    break;
                case 0:
                    freqOneButton->setToggleState(true, sendNotificationSync);
                    break;
                }
            }     
        }
    }

    void Editor::refreshLowCut()
    {
        auto p = static_cast<Node*>(getProcessor());
        lowCutEditable->setText(String(p->lowCut), dontSendNotification);
    }

    void Editor::refreshHighCut()
    {
        auto p = static_cast<Node*>(getProcessor());
        highCutEditable->setText(String(p->highCut), dontSendNotification);
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