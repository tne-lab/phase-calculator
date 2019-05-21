/*
   Previously called mempar()
   Originally in FORTRAN, hence the array offsets of 1, Yuk.
   Original code from Kay, 1988, appendix 8D.

   Perform Burg's Maximum Entropy AR parameter estimation
   outputting (or not) successive models en passant. Sourced from Alex Sergejew

   Two small changes made by NH in November 1998:
   tstarz.h no longer included, just say "typedef double REAL" instead
   Declare ar by "REAL **ar" instead of "REAL ar[MAXA][MAXA]

   Further "cleaning" by Paul Bourke.....for personal style only.

   Converted to zero-based arrays by Paul Sanders, June 2007

   Simplified and 'g' removed, plus class wrapper added by Ethan Blackwood, 2018
*/

#ifndef AR_MODELER_H_INCLUDED
#define AR_MODELER_H_INCLUDED

#include <BasicJuceHeader.h>

namespace PhaseCalculator
{
    class ARModeler {
    public:
        ARModeler(int order = 1, int length = 2, int strideIn = 1, bool* success = nullptr)
        {
            bool s = setParams(order, length, strideIn);
            if (success != nullptr)
            {
                *success = s;
            }
        }

        ~ARModeler() { }

        void reset()
        {
            hasBeenUsed = false;
        }

        bool hasBeenFit()
        {
            return hasBeenUsed;
        }

        // returns true if successful.
        bool setParams(int order, int length, int strideIn)
        {
            int newStridedLength = calcStridedLength(length, strideIn);
            if (order < 1 || newStridedLength < order + 1)
            {
                jassertfalse;
                return false;
            }
            arOrder = order;
            inputLength = length;
            stride = strideIn;
            stridedLength = newStridedLength;
            reallocateStorage();
            return true;
        }

        void getModel(Array<double, CriticalSection>& modelOut)
        {
            jassert(hasBeenUsed);

            // copies using internal mutex
            modelOut = coefficients;
        }

        void fitModel(const Array<double>& j_inputseries_reverse)
        {
            jassert(j_inputseries_reverse.size() >= inputLength);
            double t1, t2;
            int n;

            // get raw pointers to improve performance
            const double* inputseries_last = j_inputseries_reverse.begin() + inputLength - 1;
            double* coef = j_coef_temp.begin();
            double* per = j_per.begin();
            double* pef = j_pef.begin();
            double* h = j_h.begin();

            // reset per and pef
            resetPredictionError();

            for (n = 1; n <= arOrder; n++)
            {
                double sn = 0.0;
                double sd = 0.0;
                int j;
                int jj = stridedLength - n;

                for (j = 0; j < jj; j++)
                {
                    t1 = inputseries_last[-stride * (j + n)] + pef[j];
                    t2 = inputseries_last[-stride * j] + per[j];
                    sn -= 2.0 * t1 * t2;
                    sd += (t1 * t1) + (t2 * t2);
                }

                t1 = sn / sd;
                coef[n - 1] = t1;
                if (n != 1)
                {
                    for (j = 1; j < n; j++)
                        h[j - 1] = coef[j - 1] + t1 * coef[n - j - 1];
                    for (j = 1; j < n; j++)
                        coef[j - 1] = h[j - 1];
                    jj--;
                }

                for (j = 0; j < jj; j++)
                {
                    per[j] = per[j] + t1 * pef[j] + t1 * inputseries_last[-stride * (j + n)];
                    pef[j] = pef[j + 1] + t1 * per[j + 1] + t1 * inputseries_last[-stride * (j + 1)];
                }
            }

            // save (using built-in mutex)
            coefficients.swapWith(j_coef_temp);

            hasBeenUsed = true;
        }

    private:

        void reallocateStorage()
        {
            j_h.resize(arOrder - 1);
            j_per.resize(stridedLength);
            j_pef.resize(stridedLength);
            j_coef_temp.resize(arOrder);
            coefficients.resize(arOrder);
            resetPredictionError();
        }

        void resetPredictionError()
        {
            FloatVectorOperations::clear(j_per.begin(), stridedLength);
            FloatVectorOperations::clear(j_pef.begin(), stridedLength);
        }

        static int calcStridedLength(int inputLength, int stride)
        {
            jassert(stride > 0);
            return (inputLength + (stride - 1)) / stride;
        }

        int arOrder;
        int inputLength;
        int stridedLength;
        int stride;
        Array<double> j_per;
        Array<double> j_pef;
        Array<double> j_h;

        Array<double, CriticalSection> j_coef_temp;
        Array<double, CriticalSection> coefficients;

        bool hasBeenUsed = false;

        JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(ARModeler);
    };
}

#endif // AR_MODELER_H_INCLUDED
