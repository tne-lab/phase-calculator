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

#include "../../../JuceLibraryCode/JuceHeader.h"

class ARModeler {
public:
    ARModeler() = delete;
    ARModeler(int order, int length, bool* success = nullptr) : arOrder(order), inputLength(length)
    {
        if (order < 1 || order >= length)
        {
            // invalid length/order combination
            jassertfalse;
            if (success != nullptr) { *success = false; }
            arOrder = 1;
            inputLength = 2;
        }
        else if (success != nullptr) { *success = true; }
        reallocateStorage();
    }

    ~ARModeler() { }

    // returns true if successful.
    bool setOrder(int order)
    {
        if (order < 1 || order >= inputLength)
        {
            jassertfalse;
            return false;
        }
        arOrder = order;
        reallocateStorage();
        return true;
    }

    // returns true if successful.
    bool setInputLength(int length)
    {
        if (length <= arOrder)
        {
            jassertfalse;
            return false;
        }
        inputLength = length;
        reallocateStorage();
        return true;
    }

    void fitModel(const Array<double>& inputseries, Array<double>& coef)
    {
        jassert(inputseries.size() == inputLength);
        jassert(coef.size() == arOrder);
        double t1, t2;
        int n;

        // reset per and pef
        resetPredictionError();

        for (n = 1; n <= arOrder; n++)
        {
            double sn = 0.0;
            double sd = 0.0;
            int j;
            int jj = inputLength - n;

            for (j = 0; j < jj; j++)
            {
                t1 = inputseries[j + n] + pef[j];
                t2 = inputseries[j] + per[j];
                sn -= 2.0 * t1 * t2;
                sd += (t1 * t1) + (t2 * t2);
            }

            t1 = sn / sd;
            coef.setUnchecked(n - 1, t1);
            if (n != 1)
            {
                for (j = 1; j < n; j++)
                    h.setUnchecked(j - 1, coef[j - 1] + t1 * coef[n - j - 1]);
                for (j = 1; j < n; j++)
                    coef.setUnchecked(j - 1, h[j - 1]);
                jj--;
            }

            for (j = 0; j < jj; j++)
            {
                per.setUnchecked(j, per[j] + t1 * pef[j] + t1 * inputseries[j + n]);
                pef.setUnchecked(j, pef[j + 1] + t1 * per[j + 1] + t1 * inputseries[j + 1]);
            }
        }
    }

private:

    void reallocateStorage()
    {
        h.resize(arOrder - 1);
        resetPredictionError();
    }

    void resetPredictionError()
    {
        per.clearQuick();
        per.insertMultiple(0, 0, inputLength);
        pef.clearQuick();
        pef.insertMultiple(0, 0, inputLength);
    }

    int arOrder;
    int inputLength;
    Array<double> per;
    Array<double> pef;
    Array<double> h;
};

#endif AR_MODELER_H_INCLUDED