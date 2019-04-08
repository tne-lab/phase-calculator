/*
------------------------------------------------------------------

This file is part of a plugin for the Open Ephys GUI
Copyright (C) 2018 Translational NeuroEngineering Lab

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

#ifndef H_TRANSFORMERS_H_INCLUDED
#define H_TRANSFORMERS_H_INCLUDED

#include <map>

#include "../../../JuceLibraryCode/JuceHeader.h"

/*

Defines the Hilbert transformers appropriate to use for each frequency band.
(The actual coeficcients and other values are in the corresponding cpp file.)
- BAND_NAME:    display name for each frequency band.
- VALID_BAND:   range of frequencies appropriate to use with each transformer
- DEFAULT_BAND: band filled in by default when selecting each transformer
- EXTREMA:      locations of local extrema of the magnitude response within VALID_BAND (if any).
-               used to find the maximum and minimum response given a band of interest.
- DELAY:        group delay of the transformer = # of samples to predict
                also, the number of unique nonzero coefficients.
- TRANSFORMER:  first DELAY coefficients for each filter.
                the remaining DELAY+1 are 0, followed by the leading coefficients
                again, negated and in reverse order.
*/

enum Band
{
    ALPHA_THETA = 0,
    BETA,
    LOW_GAM,
    MID_GAM,
    HIGH_GAM,
    NUM_BANDS
};

namespace Hilbert
{
    const int FS = 500;

    extern const std::map<int, String> BAND_NAME;

    // each is a pair (lower limit, upper limit)
    extern const std::map<int, Array<float>> VALID_BAND;

    // each is a pair (low cut, high cut)
    extern const std::map<int, Array<float>> DEFAULT_BAND;

    extern const std::map<int, Array<float>> EXTREMA;

    // samples of group delay (= order of filter / 2)
    extern const std::map<int, int> DELAY;
    
    // contain the first delay[band] coefficients; the rest are redundant and can be inferred
    extern const std::map<int, const double*> TRANSFORMER;
}

#endif // H_TRANSFORMERS_H_INCLUDED