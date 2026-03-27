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

namespace PhaseCalculator
{
    Editor::Editor(Node* parentNode)
        : VisualizerEditor(parentNode, "Event Phase Plot", 300)
    {
        int filterWidth = 105;

        // make the canvas now, so that restoring its parameters always works.
        canvas = std::make_unique<Canvas>(parentNode);

        addComboBoxParameterEditor(Parameter::STREAM_SCOPE, "freq_range", 10, 25);

        addTextBoxParameterEditor(Parameter::STREAM_SCOPE, "low_cut", filterWidth, 25);

        addTextBoxParameterEditor(Parameter::STREAM_SCOPE, "high_cut", filterWidth, 65);

        addTextBoxParameterEditor(Parameter::STREAM_SCOPE, "ar_refresh", filterWidth + 90, 25);

        addTextBoxParameterEditor(Parameter::STREAM_SCOPE, "ar_order", filterWidth + 90, 65);

        addSelectedChannelsParameterEditor(Parameter::STREAM_SCOPE, "Channels", 10, 65);

        addComboBoxParameterEditor(Parameter::STREAM_SCOPE, "output", 10, 88);
    
        addCustomParameterEditor( new PhaseCalculatorPathParameterEditor(parentNode->getParameter("filter_config")), filterWidth, 108);

       for (auto ed : parameterEditors)
         {
          if( ed->getParameterName() == "filter_config")
              {ed->setLayout (ParameterEditor::Layout::nameHidden);
               ed->setSize(180,20); }
          else
              {ed->setLayout (ParameterEditor::Layout::nameOnTop);
               ed->setSize (90, 40); }
         }
    }

    Editor::~Editor() {}


    Visualizer* Editor::createNewCanvas()
    {
        return canvas.get();
    }

    void Editor::selectedStreamHasChanged()
    {
        LOGD("[PhaseCalc] Selected stream has changed to: ", getCurrentStream());
        Node* processor = (Node*)getProcessor();
        processor->setSelectedStream(getCurrentStream());

        // inform the canvas about selected stream updates
        updateVisualizer();
    }

    void ClearButton::paintButton (Graphics& g, bool isMouseOverButton, bool isButtonDown)
    {
    g.fillAll (findColour (ThemeColours::widgetBackground).contrasting (0.05f));

    int xoffset = (getWidth() - 14) / 2;
    int yoffset = (getHeight() - 14) / 2;

    if (isMouseOverButton)
    {
        g.setColour (findColour (ThemeColours::defaultFill).withAlpha (0.5f));
        g.fillRoundedRectangle (xoffset, yoffset, 14, 14, 4.0f);
    }

    g.setColour (findColour (ThemeColours::defaultText));

    Path path;
    path.addLineSegment (Line<float> (2, 2, 8, 8), 0.0f);
    path.addLineSegment (Line<float> (2, 8, 8, 2), 0.0f);

    path.applyTransform (AffineTransform::translation (xoffset + 2, yoffset + 2));

    g.strokePath (path, PathStrokeType (1.0f));
    }

    PhaseCalculatorPathParameterEditor::PhaseCalculatorPathParameterEditor (Parameter* param, int rowHeightPixels, int rowWidthPixels) : ParameterEditor (param)
    {
      jassert (param->getType() == Parameter::PATH_PARAM);

      setBounds (0, 0, rowWidthPixels, rowHeightPixels);

      button = std::make_unique<TextButton> ("Browse");
      button->setName (param->getKey());
      button->addListener (this);
      button->setClickingTogglesState (false);
      button->setTooltip (param->getValueAsString());
      button->addMouseListener (this, true);
      addAndMakeVisible (button.get());

      label = std::make_unique<Label> ("Parameter name", param->getDisplayName()); 
      Font labelFont = FontOptions ("Inter", "Regular", int (0.75 * rowHeightPixels));
      label->setFont (labelFont);
      label->setJustificationType (Justification::left);
      addAndMakeVisible (label.get());

      clearButton = std::make_unique<ClearButton>();
      clearButton->addListener (this);
      clearButton->addMouseListener (this, true);
      clearButton->setTooltip ("Set to default");
      addChildComponent (clearButton.get());

      int width = rowWidthPixels;

      label->setBounds (width / 2 + 4, 0, width / 2, rowHeightPixels);
      button->setBounds (0, 0, width / 2, rowHeightPixels);
      clearButton->setBounds ((width) / 2 - rowHeightPixels - 1, 1, rowHeightPixels - 2, rowHeightPixels - 2);

      editor = (Component*) button.get();
     }

    void PhaseCalculatorPathParameterEditor::buttonClicked (Button* button_)
{
    if (button_ == button.get())
    {
       auto* pathParam = static_cast<PathParameter*>(param);
        
        String dialogBoxTitle = param->getDescription().isEmpty() ?  "Select a JSON file" : param->getDescription();
        String validFilePatterns = "*." + pathParam->getValidFilePatterns().joinIntoString(";*.");

        FileChooser chooser (dialogBoxTitle, File(), validFilePatterns);

        if (chooser.browseForFileToOpen())
        {
            param->setNextValue (chooser.getResult().getFullPathName());
            updateView();
        }
    }
    else if (button_ == clearButton.get())
    {
        param->setNextValue ("None");
        clearButton->setVisible (false);
        updateView();
    }
}

void PhaseCalculatorPathParameterEditor::updateView()
{
    if (param == nullptr)
        button->setEnabled (false);
    else
        button->setEnabled (true);

    if (param)
    {
        String value = param->getValueAsString();
        button->setButtonText (value);
        if (! ((PathParameter*) param)->isValid())
        {
            button->setColour (TextButton::textColourOnId, Colours::red);
            button->setColour (TextButton::textColourOffId, Colours::red);
        }
        else
        {
            // Remove the red colour if it was previously set
            button->removeColour (TextButton::textColourOnId);
            button->removeColour (TextButton::textColourOffId);
        }
        if (value == "None")
        {
            button->setButtonText ("default");
            button->setTooltip ("Select a JSON file to override Hilbert coefficients for the available frequency bands");
        }
        else
        {
            button->setTooltip (value);
        }
    }
}

void PhaseCalculatorPathParameterEditor::resized()
{
    updateBounds();

    if (button != nullptr)
    {
        clearButton->setBounds (button->getRight() - button->getHeight(),
                                button->getY() + 1,
                                button->getHeight() - 2,
                                button->getHeight() - 2);
    }
}

void PhaseCalculatorPathParameterEditor::mouseEnter (const MouseEvent& event)
{
    if (param->getValueAsString() != "None" && button->isEnabled())
    {
        clearButton->setVisible (true);
    }
}

void PhaseCalculatorPathParameterEditor::mouseExit (const MouseEvent& event)
{
    clearButton->setVisible (false);
}
}