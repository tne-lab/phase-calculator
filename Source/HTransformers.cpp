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

#include "HTransformers.h"

namespace PhaseCalculator
{
    namespace
    {
        /** Helper to allow assigning to each array one band at a time (within the constructor) **/
        struct HilbertInfo
        {
            const String gamma{ L"\u03b3" };
            const String beta{ L"\u03b2" };
            const String alphaTheta{ L"\u03b1/\u03b8" };

            String bandName[NUM_BANDS];
            Array<float> validBand[NUM_BANDS];
            Array<float> defaultBand[NUM_BANDS];
            Array<float> extrema[NUM_BANDS];
            int delay[NUM_BANDS];
            Array<double> transformer[NUM_BANDS];

            static String validBandToString(const Array<float>& band)
            {
                jassert(band.size() == 2);
                return " (" + String(band[0]) + "-" + String(band[1]) + " Hz)";
            }

            HilbertInfo()
            {
                validBand[ALPHA_THETA] = Array<float>({ 4, 18 });
                bandName[ALPHA_THETA] = alphaTheta + validBandToString(validBand[ALPHA_THETA]);
                defaultBand[ALPHA_THETA] = Array<float>({ 4, 8 });
                extrema[ALPHA_THETA] = Array<float>(/* none */);
                delay[ALPHA_THETA] = 9;
                // from Matlab: firpm(18, [4 246]/250, [1 1], 'hilbert')
                transformer[ALPHA_THETA] = Array<double>({
                    -0.28757250783614413,
                    0.000027647225074994485,
                    -0.094611325643268351,
                    -0.00025887439499763831,
                    -0.129436276914844,
                    -0.0001608427426424053,
                    -0.21315096860055227,
                    -0.00055322197399797961,
                    -0.63685698210351149
                });


                validBand[BETA] = Array<float>({ 10, 40 });
                bandName[BETA] = beta + validBandToString(validBand[BETA]);
                defaultBand[BETA] = Array<float>({ 12, 30 });
                extrema[BETA] = Array<float>({ 21.5848 });
                delay[BETA] = 9;
                // from Matlab: firpm(18, [12 30 40 240]/250, [1 1 0.7 0.7], 'hilbert')
                transformer[BETA] = Array<double>({
                    -0.099949575596234311,
                    -0.020761484963254036,
                    -0.080803573080958854,
                    -0.027365064225587619,
                    -0.11114477443975329,
                    -0.025834076852645271,
                    -0.16664116044989324,
                    -0.015661948619847599,
                    -0.45268524264113719
                });


                validBand[LOW_GAM] = Array<float>({ 30, 55 });
                bandName[LOW_GAM] = "Lo " + gamma + validBandToString(validBand[LOW_GAM]);
                defaultBand[LOW_GAM] = Array<float>({ 30, 55 });
                extrema[LOW_GAM] = Array<float>({ 43.3609 });
                delay[LOW_GAM] = 2;
                // from Matlab: firls(4, [30 55]/250, [1 1], 'hilbert')
                transformer[LOW_GAM] = Array<double>({
                    -1.5933788446351915,
                    1.7241339075391682
                });


                validBand[MID_GAM] = Array<float>({ 40, 90 });
                bandName[MID_GAM] = "Mid " + gamma + validBandToString(validBand[MID_GAM]);
                defaultBand[MID_GAM] = Array<float>({ 40, 90 });
                extrema[MID_GAM] = Array<float>({ 64.4559 });
                delay[MID_GAM] = 2;
                // from Matlab: firls(4, [35 90]/250, [1 1], 'hilbert')
                transformer[MID_GAM] = Array<double>({
                    -0.487176162115735,
                    -0.069437334858668653
                });


                validBand[HIGH_GAM] = Array<float>({ 60, 200 });
                bandName[HIGH_GAM] = "Hi " + gamma + validBandToString(validBand[HIGH_GAM]);
                defaultBand[HIGH_GAM] = Array<float>({ 70, 150 });
                extrema[HIGH_GAM] = Array<float>({ 81.6443, 123.1104, 169.3574 });
                delay[HIGH_GAM] = 3;
                // from Matlab: firls(6, [60 200]/250, [1 1], 'hilbert')
                transformer[HIGH_GAM] = Array<double>({
                    -0.10383410506573287,
                    0.0040553935691102303,
                    -0.59258484603659545
                });
            }
        };

        static const HilbertInfo hilbertInfo; // instantiates all the data through the constructor
    }

    namespace Hilbert
    {
        // exported constants
        extern const String* const bandName = hilbertInfo.bandName;

        extern const Array<float>* const validBand = hilbertInfo.validBand;

        extern const Array<float>* const defaultBand = hilbertInfo.defaultBand;

        extern const Array<float>* const extrema = hilbertInfo.extrema;

        extern const int* const delay = hilbertInfo.delay;

        extern const Array<double>* const transformer = hilbertInfo.transformer;
    }
}