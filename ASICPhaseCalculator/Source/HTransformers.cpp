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
			const String delta{ L"\u03b4" };
			const String gamma{ L"\u03b3" };
			const String beta{ L"\u03b2" };
			const String Theta{ L"\u03b8" };
			const String alpha{ L"\u03b1" };

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
				// Added by sumedh and modified
				int del = 7;
				Array<double> htcoeff = Array<double>({ -0.433593750000000,0.0,-0.136718750000000,0.0,-0.216796875000000,0.0,-0.638671875000000 });

				/*delta band*/
				validBand[DELTA] = Array<float>({ 1, 4 });
				bandName[DELTA] = delta + validBandToString(validBand[DELTA]);
				defaultBand[DELTA] = Array<float>({ 1, 4 });
				extrema[DELTA] = Array<float>(/* none */);
				delay[DELTA] = del;
				transformer[DELTA] = htcoeff;

				/*theta band*/
				validBand[THETA] = Array<float>({ 4, 8 });
				bandName[THETA] = Theta + validBandToString(validBand[THETA]);
				defaultBand[THETA] = Array<float>({ 4, 8 });
				extrema[THETA] = Array<float>(/* none */);
				delay[THETA] = del;
				transformer[THETA] = htcoeff;

				/*alpha band*/
				validBand[ALPHA] = Array<float>({ 8, 13 });
				bandName[ALPHA] = alpha + validBandToString(validBand[ALPHA]);
				defaultBand[ALPHA] = Array<float>({ 8, 13 });
				extrema[ALPHA] = Array<float>(/* none */);
				delay[ALPHA] = del;
				transformer[ALPHA] = htcoeff;

				/*beta band*/
				validBand[BETA] = Array<float>({ 12, 30 });
				bandName[BETA] = beta + validBandToString(validBand[BETA]);
				defaultBand[BETA] = Array<float>({ 12, 30 });
				extrema[BETA] = Array<float>(/* none */);
				delay[BETA] = del;
				transformer[BETA] = htcoeff;

				/*lower gamma band*/
				validBand[LOW_GAM] = Array<float>({ 30, 55 });
				bandName[LOW_GAM] = "Lo " + gamma + validBandToString(validBand[LOW_GAM]);
				defaultBand[LOW_GAM] = Array<float>({ 30, 55 });
				extrema[LOW_GAM] = Array<float>(/* none */);
				delay[LOW_GAM] = del;
				transformer[LOW_GAM] = htcoeff;

				/*mid gamma band*/
				validBand[MID_GAM] = Array<float>({ 40, 90 });
				bandName[MID_GAM] = "Mid " + gamma + validBandToString(validBand[MID_GAM]);
				//defaultBand[MID_GAM] = Array<float>({ 40, 90 });
				extrema[MID_GAM] = Array<float>(/* none */);
				extrema[MID_GAM] = Array<float>({ 64.4559 });
				delay[MID_GAM] = del;
				transformer[MID_GAM] = htcoeff;

				/*high gamma band*/
				validBand[HIGH_GAM] = Array<float>({ 60, 200 });
				bandName[HIGH_GAM] = "Hi " + gamma + validBandToString(validBand[HIGH_GAM]);
				defaultBand[HIGH_GAM] = Array<float>({ 60, 200 });
				extrema[HIGH_GAM] = Array<float>(/* none */);
				//extrema[HIGH_GAM] = Array<float>({ 81.6443, 123.1104, 169.3574 });
				delay[HIGH_GAM] = del;
				transformer[HIGH_GAM] = htcoeff;

				/*Ripple band*/
				validBand[RIPPLE] = Array<float>({ 150, 250 });
				bandName[RIPPLE] = "RIP " + validBandToString(validBand[RIPPLE]);
				defaultBand[RIPPLE] = Array<float>({ 150, 250 });
				extrema[RIPPLE] = Array<float>(/* none */);
				delay[RIPPLE] = del;
				transformer[RIPPLE] = htcoeff;
			}
		};

		static const HilbertInfo hilbertInfo; // instantiates all the data through the constructor

		struct BandpassfiltInfo
		{
			String bandName[NUM_BANDS];
			int delay[NUM_BANDS];
			Array<double> transformer[NUM_BANDS];
			BandpassfiltInfo()
			{
				int del = 42;
				/*delta band*/
				delay[DELTA] = del;
				transformer[DELTA] = Array<double>({
				0.0761718750000000,0.0800781250000000,0.0957031250000000,0.123046875000000,0.158203125000000,0.201171875000000,0.253906250000000,0.312500000000000,0.376953125000000,0.443359375000000,0.513671875000000,0.585937500000000,0.654296875000000,0.722656250000000,0.785156250000000,0.841796875000000,0.892578125000000,0.933593750000000,0.964843750000000,0.986328125000000,0.998046875000000,0.998046875000000,0.986328125000000,0.964843750000000,0.933593750000000,0.892578125000000,0.841796875000000,0.785156250000000,0.722656250000000,0.654296875000000,0.585937500000000,0.513671875000000,0.443359375000000,0.376953125000000,0.312500000000000,0.253906250000000,0.201171875000000,0.158203125000000,0.123046875000000,0.0957031250000000,0.0800781250000000,0.0761718750000000 });

				/*theta band*/
				delay[THETA] = del;
				transformer[THETA] = Array<double>({
				0.0566406250000000,0.0625000000000000,0.0761718750000000,0.0996093750000000,0.132812500000000,0.171875000000000,0.222656250000000,0.277343750000000,0.341796875000000,0.408203125000000,0.480468750000000,0.552734375000000,0.626953125000000,0.697265625000000,0.765625000000000,0.826171875000000,0.880859375000000,0.925781250000000,0.960937500000000,0.986328125000000,0.998046875000000,0.998046875000000,0.986328125000000,0.960937500000000,0.925781250000000,0.880859375000000,0.826171875000000,0.765625000000000,0.697265625000000,0.626953125000000,0.552734375000000,0.480468750000000,0.408203125000000,0.341796875000000,0.277343750000000,0.222656250000000,0.171875000000000,0.132812500000000,0.0996093750000000,0.0761718750000000,0.0625000000000000,0.0566406250000000 });

				/*alpha band*/
				delay[ALPHA] = del;
				transformer[ALPHA] = Array<double>({
				 0.0175781250000000,0.0234375000000000,0.0351562500000000,0.0507812500000000,0.0742187500000000,0.107421875000000,0.148437500000000,0.199218750000000,0.259765625000000,0.326171875000000,0.400390625000000,0.478515625000000,0.558593750000000,0.638671875000000,0.716796875000000,0.789062500000000,0.855468750000000,0.910156250000000,0.953125000000000,0.982421875000000,0.998046875000000,0.998046875000000,0.982421875000000,0.953125000000000,0.910156250000000,0.855468750000000,0.789062500000000,0.716796875000000,0.638671875000000,0.558593750000000,0.478515625000000,0.400390625000000,0.326171875000000,0.259765625000000,0.199218750000000,0.148437500000000,0.107421875000000,0.0742187500000000,0.0507812500000000,0.0351562500000000,0.0234375000000000,0.0175781250000000 });

				/*beta band*/
				delay[BETA] = del;
				transformer[BETA] = Array<double>({
				-0.0566406250000000,-0.0585937500000000,-0.0644531250000000,-0.0722656250000000,-0.0800781250000000,-0.0839843750000000,-0.0781250000000000,-0.0605468750000000,-0.0273437500000000,0.0234375000000000,0.0917968750000000,0.175781250000000,0.277343750000000,0.388671875000000,0.505859375000000,0.623046875000000,0.734375000000000,0.833984375000000,0.912109375000000,0.968750000000000,0.998046875000000,0.998046875000000,0.968750000000000,0.912109375000000,0.833984375000000,0.734375000000000,0.623046875000000,0.505859375000000,0.388671875000000,0.277343750000000,0.175781250000000,0.0917968750000000,0.0234375000000000,-0.0273437500000000,-0.0605468750000000,-0.0781250000000000,-0.0839843750000000,-0.0800781250000000,-0.0722656250000000,-0.0644531250000000,-0.0585937500000000,-0.0566406250000000 });

				/*low gamma band*/
				delay[LOW_GAM] = del;
				transformer[LOW_GAM] = Array<double>({
					0.0351562500000000,0.0273437500000000,0.0156250000000000,-0.00390625000000000,-0.0371093750000000,-0.0878906250000000,-0.156250000000000,-0.238281250000000,-0.322265625000000,-0.396484375000000,-0.443359375000000,-0.447265625000000,-0.398437500000000,-0.289062500000000,-0.125000000000000,0.0839843750000000,0.318359375000000,0.552734375000000,0.761718750000000,0.916015625000000,0.998046875000000,0.998046875000000,0.916015625000000,0.761718750000000,0.552734375000000,0.318359375000000,0.0839843750000000,-0.125000000000000,-0.289062500000000,-0.398437500000000,-0.447265625000000,-0.443359375000000,-0.396484375000000,-0.322265625000000,-0.238281250000000,-0.156250000000000,-0.0878906250000000,-0.0371093750000000,-0.00390625000000000,0.0156250000000000,0.0273437500000000,0.0351562500000000 });

				/*mid gamma band*/
				delay[MID_GAM] = del;
				transformer[MID_GAM] = Array<double>({
					0.00195312500000000,0.0,0.00195312500000000,0.0117187500000000,0.0312500000000000,0.0566406250000000,0.0839843750000000,0.0937500000000000,0.0703125000000000,-0.00390625000000000,-0.132812500000000,-0.298828125000000,-0.466796875000000,-0.582031250000000,-0.595703125000000,-0.476562500000000,-0.220703125000000,0.128906250000000,0.501953125000000,0.818359375000000,0.998046875000000,0.998046875000000,0.818359375000000,0.501953125000000,0.128906250000000,-0.220703125000000,-0.476562500000000,-0.595703125000000,-0.582031250000000,-0.466796875000000,-0.298828125000000,-0.132812500000000,-0.00390625000000000,0.0703125000000000,0.0937500000000000,0.0839843750000000,0.0566406250000000,0.0312500000000000,0.0117187500000000,0.00195312500000000,0.0,0.00195312500000000 });

				/*High gamma band*/
				delay[HIGH_GAM] = del;
				transformer[HIGH_GAM] = Array<double>({
					-0.00195312500000000,-0.00781250000000000,-0.0117187500000000,-0.00195312500000000,0.0117187500000000,0.0175781250000000,0.00390625000000000,0.0,0.0390625000000000,0.0917968750000000,0.0820312500000000,-0.0117187500000000,-0.0859375000000000,-0.0371093750000000,0.0468750000000000,-0.0546875000000000,-0.392578125000000,-0.640625000000000,-0.390625000000000,0.341796875000000,0.998046875000000,0.998046875000000,0.341796875000000,-0.390625000000000,-0.640625000000000,-0.392578125000000,-0.0546875000000000,0.0468750000000000,-0.0371093750000000,-0.0859375000000000,-0.0117187500000000,0.0820312500000000,0.0917968750000000,0.0390625000000000,0.0,0.00390625000000000,0.0175781250000000,0.0117187500000000,-0.00195312500000000,-0.0117187500000000,-0.00781250000000000,-0.00195312500000000 });

				/*Ripple band*/
				delay[RIPPLE] = del;
				transformer[RIPPLE] = Array<double>({
					 0.00195312500000000,-0.00195312500000000,0.00195312500000000,0.0195312500000000,0.00976562500000000,-0.0390625000000000,-0.0527343750000000,0.0234375000000000,0.0800781250000000,0.0195312500000000,-0.0234375000000000,0.0292968750000000,-0.0410156250000000,-0.251953125000000,-0.121093750000000,0.449218750000000,0.580078125000000,-0.269531250000000,-0.998046875000000,-0.335937500000000,0.921875000000000,0.921875000000000,-0.335937500000000,-0.998046875000000,-0.269531250000000,0.580078125000000,0.449218750000000,-0.121093750000000,-0.251953125000000,-0.0410156250000000,0.0292968750000000,-0.0234375000000000,0.0195312500000000,0.0800781250000000,0.0234375000000000,-0.0527343750000000,-0.0390625000000000,0.00976562500000000,0.0195312500000000,0.00195312500000000,-0.00195312500000000,0.00195312500000000 });

			}
		};

		static const BandpassfiltInfo bandpassfiltInfo; // instantiates all the data through the constructor

		struct LowpassfiltInfo
		{
			String bandName[NUM_BANDS];
			int delay[NUM_BANDS];
			Array<double> transformer[NUM_BANDS];

			LowpassfiltInfo()
			{
				int del = 25;
				Array<double> lpfcoeff = Array<double>({ 0.0,0.00195312500000000,0.00585937500000000,0.00781250000000000,0.0,-0.0332031250000000,-0.0800781250000000,-0.0917968750000000,0.0,0.238281250000000,0.574218750000000,0.875000000000000,0.998046875000000,0.875000000000000,0.574218750000000,0.238281250000000,0.0,-0.0917968750000000,-0.0800781250000000,-0.0332031250000000,0.0,0.00781250000000000,0.00585937500000000,0.00195312500000000,0.0 });
				/*delta band*/
				delay[DELTA] = del;
				transformer[DELTA] = lpfcoeff;
				/*theta band*/
				delay[THETA] = del;
				transformer[THETA] = lpfcoeff;
				/*alpha band*/
				delay[ALPHA] = del;
				transformer[ALPHA] = lpfcoeff;
				/*beta band*/
				delay[BETA] = del;
				transformer[BETA] = lpfcoeff;
				/*low gamma band*/
				delay[LOW_GAM] = del;
				transformer[LOW_GAM] = lpfcoeff;
				/*mid gamma band*/
				delay[MID_GAM] = del;
				transformer[MID_GAM] = lpfcoeff;
				/*High gamma band*/
				delay[HIGH_GAM] = del;
				transformer[HIGH_GAM] = lpfcoeff;
				/*Ripple band*/
				delay[RIPPLE] = del;
				transformer[RIPPLE] = lpfcoeff;
			}
		};

		static const LowpassfiltInfo lowpassfiltInfo; // instantiates all the data through the constructor
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

	namespace bandpassfilt
	{
		// exported constants
		extern const String* const bandName = bandpassfiltInfo.bandName;

		extern const int* const delay = bandpassfiltInfo.delay;

		extern const Array<double>* const transformer = bandpassfiltInfo.transformer;
	}
	namespace lowpassfilt
	{
		// exported constants
		extern const String* const bandName = lowpassfiltInfo.bandName;

		extern const int* const delay = lowpassfiltInfo.delay;

		extern const Array<double>* const transformer = lowpassfiltInfo.transformer;
	}
}