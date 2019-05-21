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

#include <JuceHeader.h> // assertions, etc.
#include <complex>
#include <algorithm>    // reverse array
#include <utility>

#include "fftw3.h" // Fast Fourier Transform library

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

    FFTWArray(const FFTWArray& other)
        : FFTWArray(other.length)
    {
        // delegate to copy assignment
        *this = other;
    }

    // Move to new FFTWArray. Needed to store in an array
    FFTWArray(FFTWArray&& other)
    {
        // delegate to move assignment
        *this = std::move(other);
    }

    FFTWArray& operator=(const FFTWArray& other)
    {
        if (this != &other)
        {
            resize(other.length);
            copyFrom(other.data, other.length);
        }
        return *this;
    }

    FFTWArray& operator=(FFTWArray&& other)
    {
        if (this != &other)
        {
            data = other.data;
            length = other.length;
            other.data = nullptr;
            other.length = 0;
        }
        return *this;
    }

    // returns true if a resize actually occurred
    virtual bool resize(int newLength)
    {
        jassert(newLength >= 0);
        if (newLength != length)
        {
            length = newLength;
            fftw_free(data);
            data = reinterpret_cast<std::complex<double>*>(fftw_alloc_complex(newLength));
            return true;
        }
        return false;
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

    const double* getReadPointer(int index = 0) const
    {
        if (index < length * 2 && index >= 0)
        {
            return reinterpret_cast<const double*>(data)+index;
        }
        return nullptr;
    }

    int getLength() const
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

    // Does the part of the Hilbert transform in the frequency domain (see FFTWTransformableArray::hilbert)
    void freqDomainHilbert()
    {
        if (length <= 0) { return; }

        // normalize DC and Nyquist, normalize and double positive freqs, and set negative freqs to 0.
        int lastPosFreq = (length + 1) / 2 - 1;
        int firstNegFreq = length / 2 + 1;
        int numPosNegFreqDoubles = lastPosFreq * 2; // sizeof(complex<double>) = 2 * sizeof(double)
        bool hasNyquist = (length % 2 == 0);

        // normalize but don't double DC value
        data[0] /= length;

        // normalize and double positive frequencies
        FloatVectorOperations::multiply(reinterpret_cast<double*>(data + 1), 2.0 / length, numPosNegFreqDoubles);

        if (hasNyquist)
        {
            // normalize but don't double Nyquist frequency
            data[lastPosFreq + 1] /= length;
        }

        // set negative frequencies to 0
        FloatVectorOperations::clear(reinterpret_cast<double*>(data + firstNegFreq), numPosNegFreqDoubles);
    }
    
private:
    std::complex<double>* data;
    int length;

    JUCE_LEAK_DETECTOR(FFTWArray);
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


// copyable, movable array that can do an in-place transform
class FFTWTransformableArray : public FFTWArray
{
public:
    FFTWTransformableArray(int n = 0, unsigned int flags = FFTW_MEASURE)
        : FFTWArray(0)
        , flags(flags)
    {
        resize(n);
    }

    FFTWTransformableArray(const FFTWTransformableArray& other)
        : FFTWTransformableArray(other.getLength(), other.flags)
    {
        // delegate to copy assignment
        *this = other;
    }

    FFTWTransformableArray(FFTWTransformableArray&& other)
    {
        // delegate to move assignment
        *this = std::move(other);
    }

    FFTWTransformableArray& operator=(const FFTWTransformableArray& other)
    {
        if (this != &other)
        {
            if (flags != other.flags)
            {
                flags = other.flags;
                resize(0); // force plans to be remade
            }
            FFTWArray::operator=(other);
            // the plans will be copied if/when the overridden "resize" is called
        }
        return *this;
    }

    FFTWTransformableArray& operator=(FFTWTransformableArray&& other)
    {
        if (this != &other)
        {
            flags = other.flags;
            forwardPlan = other.forwardPlan;
            inversePlan = other.inversePlan;
            r2cPlan = other.r2cPlan;
            FFTWArray::operator=(std::move(other));
        }
        return *this;
    }

    bool resize(int newLength) override
    {
        if (FFTWArray::resize(newLength))
        {
            forwardPlan = newLength > 0 ? new FFTWPlan(newLength, this, FFTW_FORWARD, flags) : nullptr;
            inversePlan = newLength > 0 ? new FFTWPlan(newLength, this, FFTW_BACKWARD, flags) : nullptr;
            r2cPlan = newLength > 0 ? new FFTWPlan(newLength, this, flags) : nullptr;
            return true;
        }
        return false;
    }

    void fftComplex()
    {
        if (forwardPlan != nullptr)
        {
            forwardPlan->execute();
        }
    }

    void fftReal()
    {
        if (r2cPlan != nullptr)
        {
            r2cPlan->execute();
        }
    }

    void ifft()
    {
        if (inversePlan != nullptr)
        {
            inversePlan->execute();
        }
    }

    // Do Hilbert transform of real data by taking the fft, zeroing out negative frequencies,
    // and taking the ifft. The result is not technically the Hilbert transform but the
    // analytic signal, defined as x + H[x]*i. This is consistent with Matlab's 'hilbert'.
    void hilbert()
    {
        fftReal();
        freqDomainHilbert();
        ifft();
    }

private:
    unsigned int flags;
    ScopedPointer<FFTWPlan> forwardPlan;
    ScopedPointer<FFTWPlan> inversePlan;
    ScopedPointer<FFTWPlan> r2cPlan;

    JUCE_LEAK_DETECTOR(FFTWTransformableArray);
};


// version with different flags that's still default-constructible
// TODO think of a better name?
template<unsigned int f>
class FFTWTransformableArrayUsing : public FFTWTransformableArray 
{
public:
    FFTWTransformableArrayUsing(int n = 0)
        : FFTWTransformableArray(n, f)
    {}
};

#endif // FFTW_WRAPPER_H_INCLUDED