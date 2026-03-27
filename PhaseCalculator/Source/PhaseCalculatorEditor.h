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

#include <sstream>  // string parsing
#include <VisualizerEditorHeaders.h>

#include "PhaseCalculator.h"

namespace PhaseCalculator
{
    class Editor
        : public VisualizerEditor
    {
        friend class RosePlot;  // to access label updating method
    public:
        Editor(Node* parentNode);
        ~Editor();

        Visualizer* createNewCanvas() override;

    private:
        void selectedStreamHasChanged() override;

        JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(Editor);
    };

    class ClearButton : public juce::Button
    {
    public:
        /** Constructor */
        ClearButton() : Button ("Revert Dir") {}

        /** Renders the button */
        void paintButton (Graphics& g, bool isMouseOverButton, bool isButtonDown) override;
    };
    class PhaseCalculatorPathParameterEditor : public ParameterEditor,
                                  public Button::Listener
    {
      public:
      /** Constructor */
      PhaseCalculatorPathParameterEditor (Parameter* param, int rowHeightPixels = 18, int rowWidthPixels = 170);

      /** Destructor */
      virtual ~PhaseCalculatorPathParameterEditor() {}

      /** Displays the File Chooser*/
      void buttonClicked (Button* label) override;

      /** Must ensure that editor state matches underlying parameter */
      void updateView() override;

      /** Sets sub-component locations */
      void resized() override;

      /** Shows a clear button if non-default path is set*/
      void mouseEnter (const MouseEvent& event) override;

      /** Hides the clear button*/
      void mouseExit (const MouseEvent& event) override;

      private:
      std::unique_ptr<TextButton> button;
      std::unique_ptr<ClearButton> clearButton;
};
}

#endif // PHASE_CALCULATOR_EDITOR_H_INCLUDED
