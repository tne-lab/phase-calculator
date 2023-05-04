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
}

#endif // PHASE_CALCULATOR_EDITOR_H_INCLUDED
