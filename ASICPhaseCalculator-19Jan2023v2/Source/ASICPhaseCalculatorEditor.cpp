/*
------------------------------------------------------------------

This file is part of a plugin for the Open Ephys GUI
Copyright (C) 2022 Translational NeuroEngineering Laboratory, MGH

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


#include "ASICPhaseCalculator.h"
#include "ASICPhaseCalculatorEditor.h"
#include "HTransformers.h"


namespace ASICPhaseCalculator {

    Editor::~Editor() {
        node.release();                   // Releasing the Pointer to Node to avoid Deletion of ASCI Class Object as It will be done by its own destructor
        canvas.release();
    }


    Visualizer* Editor::createNewCanvas() {
        return canvas.get();
    }
    

    Editor::Editor( Node* parentNode)
        : VisualizerEditor(parentNode, 300)
    {
                
        desiredWidth = 300;                                           // Desired Width of The ASIC Plugin
        node.reset(parentNode);                                       // Setting Pointer to ASIC Processor Node
        chglistnr = std::make_unique<AllChangeListener>(this);        // Making a pointer to the Listener Class
        canvas = std::make_unique<Canvas>(parentNode);                // Creating Canvas Object and Setting pointer to it

        //addMaskChannelsParameterEditor("Channels", 10, 25);
        //addSelectedChannelsParameterEditor("Channels", 10, 25);


        bandLabel = new Label("bandL", "Freq range:");                // "Freq range:" Static Text Generator for ASIC Plugin
        bandLabel->setBounds(10, 46, 100, 20);                         
        bandLabel->setFont({ "Small Text", 12, Font::plain });
        bandLabel->setColour(Label::textColourId, Colours::darkgrey);
        addAndMakeVisible(bandLabel);


        bandBox = new ComboBox("bandB");                              // Dropdown Menue Settings for Frquency Band Selection
        for (int b = 0; b < NUM_BANDS; ++b) {
            bandBox->addItem(Hilbert::bandName[b], b + 1);            // Iterating throught the Band Names [emum ]
        }
        bandBox->setSelectedId(parentNode->getBand() + 1);            // Setting the Default Selected ID to 0+1 [theta band]
        bandBox->setTooltip(freqRangeTooltip);
        bandBox->setBounds(15, 66, 110, 20);
        bandBox->addListener(chglistnr.get());                        // Getting the Pointer to Listener Class Wrapped inside unique_ptr and Setting it up as Listener Class
        addAndMakeVisible(bandBox);


        outputModeLabel = new Label("outputModeL", "Output:");        // Generates "Output:" Label
        outputModeLabel->setBounds(10, 86, 70, 20);
        outputModeLabel->setFont({ "Small Text", 12, Font::plain });
        outputModeLabel->setColour(Label::textColourId, Colours::darkgrey);
        addAndMakeVisible(outputModeLabel);


        outputModeBox = new ComboBox("outputModeB");                  // Setting Up a Drop Down Menue for Choosing the Outputmode
        outputModeBox->addItem("PHASE", PH);                          
        outputModeBox->addItem("MAG", MAG);
        outputModeBox->addItem("PH+MAG", PH_AND_MAG);
        outputModeBox->addItem("IMAG", IM);
        outputModeBox->setSelectedId(parentNode->getOutputMode());    // Aquiring the Previous outputMode and Setting it 
        outputModeBox->setTooltip(outputModeTooltip);
        outputModeBox->setBounds(15, 106, 100, 19);
        outputModeBox->addListener(chglistnr.get());                  // Getting the Pointer to Listener Class Wrapped inside unique_ptr and Setting it up as Listener Class                   
        addAndMakeVisible(outputModeBox);


        downsampleIlabel = new Label("dsI", "Downsample I:");         // Generates label "Downsample I:"
        downsampleIlabel->setBounds(130, 25, 100, 20);
        downsampleIlabel->setFont({ "Small Text", 12, Font::bold });
        downsampleIlabel->setColour(Label::textColourId, Colours::darkgrey);
        addAndMakeVisible(downsampleIlabel);


        downsampleIlabelEdit = new Label("dsIedit");                  // Generates Edit Text Box for Entering DownSampling 1.
        downsampleIlabelEdit->setEditable(true);
        downsampleIlabelEdit->addListener(chglistnr.get());
        downsampleIlabelEdit->setBounds(230, 25, 50, 20);
        downsampleIlabelEdit->setText(String(parentNode->getDownSample1Rate()), dontSendNotification);
        downsampleIlabelEdit->setColour(Label::backgroundColourId, Colours::grey);
        downsampleIlabelEdit->setColour(Label::textColourId, Colours::white);
        addAndMakeVisible(downsampleIlabelEdit);


        downsampleIIlabel = new Label("dsII", "Downsample II:");                   // Set the Label Downsample II
        downsampleIIlabel->setBounds(130, 50, 100, 20);
        downsampleIIlabel->setFont({ "Small Text", 12, Font::bold });
        downsampleIIlabel->setColour(Label::textColourId, Colours::darkgrey);
        addAndMakeVisible(downsampleIIlabel);


        downsampleIIlabelEdit = new Label("dsIIedit");                             // Generates Edit Text Box for Entering DownSampling 2.
        downsampleIIlabelEdit->setEditable(true);
        downsampleIIlabelEdit->addListener(chglistnr.get());
        downsampleIIlabelEdit->setBounds(230, 50, 50, 19);
        downsampleIIlabelEdit->setText(String(parentNode->getDownSample2Rate()), dontSendNotification);
        downsampleIIlabelEdit->setColour(Label::backgroundColourId, Colours::grey);
        downsampleIIlabelEdit->setColour(Label::textColourId, Colours::white);
        addAndMakeVisible(downsampleIIlabelEdit);


        txtlabelTEXT = new Label("txtText", "Save Text:");                          // Adds Static Label Save Text
        txtlabelTEXT->setBounds(130, 75, 100, 20);
        txtlabelTEXT->setFont({ "Small Text", 12, Font::plain });
        txtlabelTEXT->setColour(Label::textColourId, Colours::darkgrey);
        addAndMakeVisible(txtlabelTEXT);


        istextsavebox = new ComboBox("SAVE TEXT");                                  // Dropdown Menue for Saving the Text or Not 
        istextsavebox->addItem("NO", NTEXT);
        istextsavebox->addItem("YES", YTEXT);
        istextsavebox->setSelectedId(parentNode->getIsText());
        istextsavebox->setTooltip(TxtTooltip);
        istextsavebox->setBounds(230, 75, 50, 20);
        istextsavebox->addListener(chglistnr.get());
        addAndMakeVisible(istextsavebox);


        txtlabelLAA = new Label("txtLAA", "LAA:");                                  // Adds Static Label "LAA:"
        txtlabelLAA->setBounds(130, 100, 100, 20);
        txtlabelLAA->setFont({ "Small Text", 12, Font::plain });
        txtlabelLAA->setColour(Label::textColourId, Colours::darkgrey);
        addAndMakeVisible(txtlabelLAA);


        isLAAbox = new ComboBox("LAA");                                             // Dropdown Menue for Setting LAA or Not
        isLAAbox->addItem("NO", NLAA);
        isLAAbox->addItem("YES", YLAA);
        isLAAbox->setSelectedId(parentNode->getIsLAA());
        isLAAbox->setTooltip(LAATooltip);
        isLAAbox->setBounds(230, 100, 50, 20);
        isLAAbox->addListener(chglistnr.get());
        addAndMakeVisible(isLAAbox);

        // new channels should be disabled by default
        channelSelector->paramButtonsToggledByDefault(false);

    }


    void Editor::channelChanged(int chan, bool newState){
        
        Array <int> selectedChannels = getActiveChannels();
        node->selActiveChannels = selectedChannels;
        node->recreateAllActiveChannels(selectedChannels, node->getInputChCount());
        updateVisualizer();
    }
    
    
    void Editor::startAcquisition()
    {
        bandBox->setEnabled(false);
        outputModeBox->setEnabled(false);
        channelSelector->inactivateButtons();
        //added by sumedh
        downsampleIlabelEdit->setEditable(false);
        downsampleIIlabelEdit->setEditable(false);
        istextsavebox->setEnabled(false);
        isLAAbox->setEnabled(false);

        node->setLockEditor(true); 
        node->startCalVisualizerThread();

    }

    void Editor::stopAcquisition()
    {
        bandBox->setEnabled(true);
        outputModeBox->setEnabled(true);
        channelSelector->activateButtons();
        //added by sumedh
        downsampleIlabelEdit->setEditable(true);
        downsampleIIlabelEdit->setEditable(true);
        istextsavebox->setEnabled(true);
        isLAAbox->setEnabled(true);

        node->setLockEditor(false);
        node->stopCalVisualizerThread();
    }
    
    //_________________________________________________________________________________________________________________
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                        ALL CHANGES LISTENER CLASS METHODS FOR ASIC PLUGIN
    //\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    

    AllChangeListener::AllChangeListener(Editor* editorNode) {
        editor.reset(editorNode);               // Set a Pointer to the editor Class
    }

    AllChangeListener::~AllChangeListener() {
        editor.release();                       // Releasing Unique Pointer to Editor (Prevents the Editor Object to be deleted),
                                                // The editor Object will be deleted by its own Destructor.
    }


    void AllChangeListener::comboBoxChanged(ComboBox* comboBoxThatHasChanged){
        if (editor->node->getLockEditorInfo()) {                                      // If the Editor is locked dont Update Anything [Set it Back]
            if (comboBoxThatHasChanged == editor->bandBox) {
                comboBoxThatHasChanged->setSelectedId(editor->node->getBand() + 1);
            }
            else if (comboBoxThatHasChanged == editor->outputModeBox) {
                comboBoxThatHasChanged->setSelectedId(editor->node->getOutputMode());
            }
            else if (comboBoxThatHasChanged == editor->istextsavebox) {
                comboBoxThatHasChanged->setSelectedId(editor->node->getIsText());
            }
            else if (comboBoxThatHasChanged == editor->isLAAbox) {
                comboBoxThatHasChanged->setSelectedId(editor->node->getIsLAA());
            }
            else {
                jassertfalse;
            }
            return;
        }

        if (comboBoxThatHasChanged == editor->bandBox) {
            editor->node->setBand(Band(comboBoxThatHasChanged->getSelectedId() - 1));
        }
        else if (comboBoxThatHasChanged == editor->outputModeBox) {
            editor->node->setOutputMode(OutputMode(comboBoxThatHasChanged->getSelectedId()));
        }
        else if (comboBoxThatHasChanged == editor->istextsavebox) {
            editor->node->setIsText(EisText(comboBoxThatHasChanged->getSelectedId()));
        }
        else if (comboBoxThatHasChanged == editor->isLAAbox) {
            editor->node->setIsLAA(EisLAA(comboBoxThatHasChanged->getSelectedId()));
        }
    }

   
    void AllChangeListener::labelTextChanged(Label* labelThatHasChanged) {
        if (editor->node->getLockEditorInfo()) {
            if (labelThatHasChanged == editor->downsampleIlabelEdit) {
                labelThatHasChanged->setText(String(editor->node->getDownSample1Rate()), dontSendNotification);
            }
            else if ((labelThatHasChanged == editor->downsampleIIlabelEdit)) {
                labelThatHasChanged->setText(String(editor->node->getDownSample2Rate()), dontSendNotification);
            }
            else {
                jassertfalse;
            }
            return;
        }

        int intInput;
        if (labelThatHasChanged == editor->downsampleIlabelEdit) {
            bool valid = editor->updateControl(labelThatHasChanged, 0, INT_MAX, editor->node->getDownSample1Rate(), intInput);   // Check if the Entered Number is Valid
            if (valid)
                editor->node->setDownSample1Rate(intInput);                                                                      // If Number is Valid the Update ASIC::dsIscale
            else
                labelThatHasChanged->setText(String(editor->node->getDownSample1Rate()) , dontSendNotification);                 // If Not Valied Revert Back to the Previous Value
        }
        else if (labelThatHasChanged == editor->downsampleIIlabelEdit) { 
            bool valid = editor->updateControl(labelThatHasChanged, 0, INT_MAX, editor->node->getDownSample2Rate(), intInput);   // Check if the Entered Number for ASIC::dsIIscale is valid
            if (valid)
                editor->node->setDownSample2Rate(intInput);                                                                      // If valid the updated ASIC::dsIIscale with this value
            else    
                labelThatHasChanged->setText(String(editor->node->getDownSample2Rate()), dontSendNotification);                  // If Not Valid Revret Back GUI with last Valid Value
        }
    }

    

    


}