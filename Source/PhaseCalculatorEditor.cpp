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

        addComboBoxParameterEditor("freq_range", 10, 30);

        addTextBoxParameterEditor("low_cut", filterWidth, 25);

        addTextBoxParameterEditor("high_cut", filterWidth, 75);

        addTextBoxParameterEditor("ar_refresh", filterWidth + 90, 25);

        addTextBoxParameterEditor("ar_order", filterWidth + 90, 75);

        addSelectedChannelsParameterEditor("Channels", 10, 90);

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
}