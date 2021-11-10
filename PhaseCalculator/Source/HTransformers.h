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

#include <BasicJuceHeader.h>

/*

Defines the Hilbert transformers appropriate to use for each frequency band.
(The actual coeficcients and other values are in the corresponding cpp file.)
- bandName:     display name for each frequency band.
- validBand:    range of frequencies appropriate to use with each transformer
- defaultBand:  band filled in by default when selecting each transformer
- extrema:      locations of local extrema of the magnitude response within VALID_BAND (if any).
-               used to find the maximum and minimum response given a band of interest.
- delay:        group delay of the transformer = # of samples to predict
                also, the number of unique nonzero coefficients.
- transformer:  first DELAY coefficients for each filter.
                the remaining DELAY+1 are 0, followed by the leading coefficients
                again, negated and in reverse order.
*/

namespace PhaseCalculator
{
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
        const int fs = 500;

        // Each pointer below points to an array of length NUM_BANDS.

        extern const String* const bandName;

        // each is a pair (lower limit, upper limit)
        extern const Array<float>* const validBand;

        // each is a pair (low cut, high cut)
        extern const Array<float>* const defaultBand;

        extern const Array<float>* const extrema;

        // samples of group delay (= order of filter / 2)
        extern const int* const delay;

        // contain the first delay[band] coefficients; the rest are redundant and can be inferred
        extern const Array<double>* const transformer;
    }
}

#endif // H_TRANSFORMERS_H_INCLUDED