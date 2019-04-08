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

/*
Object-oriented / RAII-friendly wrapper for relevant parts of the FFTW Fourier
transform library
*/

#ifndef FFTW_WRAPPER_H_INCLUDED
#define FFTW_WRAPPER_H_INCLUDED

#include <BasicJuceHeader.h> // assertions, etc.
#include <fftw3.h>           // Fast Fourier Transform library
#include <complex>
#include <algorithm>         // reverse array

namespace PhaseCalculator
{
    /*
    FFTW-friendly array that can hold complex or real doubles.
    */
    class FFTWArray
    {
    public:
        // creation / deletion
        FFTWArray(int complexLen = 0)
        {
            jassert(complexLen >= 0);
            length = complexLen;
            data = reinterpret_cast<std::complex<double>*>(fftw_alloc_complex(complexLen));
        }

        ~FFTWArray()
        {
            fftw_free(data);
        }

        void resize(int newLength)
        {
            jassert(newLength >= 0);
            if (newLength != length)
            {
                length = newLength;
                fftw_free(data);
                data = reinterpret_cast<std::complex<double>*>(fftw_alloc_complex(newLength));
            }
        }

        std::complex<double> getAsComplex(int i)
        {
            jassert(i >= 0 && i < length);
            return data[i];
        }

        double getAsReal(int i)
        {
            jassert(i >= 0 && i < length * 2);
            return reinterpret_cast<double*>(data)[i];
        }

        std::complex<double>* getComplexPointer(int index = 0)
        {
            if (index < length && index >= 0)
            {
                return data + index;
            }
            return nullptr;
        }

        double* getRealPointer(int index = 0)
        {
            if (index < length * 2 && index >= 0)
            {
                return reinterpret_cast<double*>(data)+index;
            }
            return nullptr;
        }

        int getLength()
        {
            return length;
        }

        // modification
        void set(int i, std::complex<double> val)
        {
            jassert(i >= 0 && i < length);
            data[i] = val;
        }

        void set(int i, double val)
        {
            jassert(i >= 0 && i < length * 2);
            reinterpret_cast<double*>(data)[i] = val;
        }

        // Reverses first reverseLength complex values (default = all)
        void reverseComplex(int reverseLength = -1)
        {
            reverseLength = reverseLength >= 0 ? reverseLength : length;
            std::complex<double>* first = data;
            std::complex<double>* last = data + reverseLength;
            std::reverse<std::complex<double>*>(first, last);
        }

        // Reverses first reverseLength real values (default = all)
        void reverseReal(int reverseLength = -1)
        {
            reverseLength = reverseLength >= 0 ? reverseLength : length * 2;
            double* first = reinterpret_cast<double*>(data);
            double* last = first + reverseLength;
            std::reverse<double*>(first, last);
        }

        /* Copies up to num elements starting at fromArr to the array starting at startInd.
        * Returns the number of elements actually copied.
        */
        int copyFrom(const std::complex<double>* fromArr, int num, int startInd = 0)
        {
            int numToCopy = jmin(num, length - startInd);
            std::complex<double>* wp = getComplexPointer(startInd);
            for (int i = 0; i < numToCopy; ++i)
            {
                wp[i] = fromArr[i];
            }

            return numToCopy;
        }

        int copyFrom(const double* fromArr, int num, int startInd = 0)
        {
            int numToCopy = jmin(num, 2 * length - startInd);
            double* wpReal = getRealPointer(startInd);
            for (int i = 0; i < numToCopy; ++i)
            {
                wpReal[i] = fromArr[i];
            }

            return numToCopy;
        }

    private:
        std::complex<double>* data;
        int length;

        JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(FFTWArray);
    };

    class FFTWPlan
    {
    public:
        // r2c constructor
        FFTWPlan(int n, FFTWArray* in, FFTWArray* out, unsigned int flags) : length(n)
        {
            double* ptr_in = in->getRealPointer();
            fftw_complex* ptr_out = reinterpret_cast<fftw_complex*>(out->getComplexPointer());
            plan = fftw_plan_dft_r2c_1d(n, ptr_in, ptr_out, flags);
        }

        // r2c in-place
        FFTWPlan(int n, FFTWArray* buf, unsigned int flags) : FFTWPlan(n, buf, buf, flags) {}

        // c2c constructor
        FFTWPlan(int n, FFTWArray* in, FFTWArray* out, int sign, unsigned int flags) : length(n)
        {
            fftw_complex* ptr_in = reinterpret_cast<fftw_complex*>(in->getComplexPointer());
            fftw_complex* ptr_out = reinterpret_cast<fftw_complex*>(out->getComplexPointer());
            plan = fftw_plan_dft_1d(n, ptr_in, ptr_out, sign, flags);
        }

        // c2c in-place
        FFTWPlan(int n, FFTWArray* buf, int sign, unsigned int flags) : FFTWPlan(n, buf, buf, sign, flags) {}

        ~FFTWPlan()
        {
            fftw_destroy_plan(plan);
        }

        void execute()
        {
            fftw_execute(plan);
        }

        int getLength()
        {
            return length;
        }

    private:
        fftw_plan plan;
        const int length;

        JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(FFTWPlan);
    };
}

#endif // FFTW_WRAPPER_H_INCLUDED