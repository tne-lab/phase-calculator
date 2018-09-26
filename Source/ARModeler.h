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
    ARModeler() : arOrder(1), inputLength(2)
    {
        reallocateStorage();
    }

    ARModeler(int order, int length, bool* success = nullptr) : ARModeler()
    {
        bool s = setParams(order, length);
        if (success != nullptr)
        {
            *success = s;
        }
    }

    ~ARModeler() { }

    // returns true if successful.
    bool setParams(int order, int length)
    {
        if (order < 1 || order >= length)
        {
            jassertfalse;
            return false;
        }
        arOrder = order;
        inputLength = length;
        reallocateStorage();
        return true;
    }

    void fitModel(const Array<double>& j_inputseries, Array<double>& j_coef)
    {
        jassert(j_inputseries.size() == inputLength);
        jassert(j_coef.size() == arOrder);
        double t1, t2;
        int n;

        // get raw pointers to improve performance
        const double* inputseries = j_inputseries.begin();
        double* coef = j_coef.begin();
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
            int jj = inputLength - n;

            for (j = 0; j < jj; j++)
            {
                t1 = inputseries[j + n] + pef[j];
                t2 = inputseries[j] + per[j];
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
                per[j] = per[j] + t1 * pef[j] + t1 * inputseries[j + n];
                pef[j] = pef[j + 1] + t1 * per[j + 1] + t1 * inputseries[j + 1];
            }
        }
    }

private:

    void reallocateStorage()
    {
        j_h.resize(arOrder - 1);
        resetPredictionError();
    }

    void resetPredictionError()
    {
        j_per.clearQuick();
        j_per.insertMultiple(0, 0, inputLength);
        j_pef.clearQuick();
        j_pef.insertMultiple(0, 0, inputLength);
    }

    int arOrder;
    int inputLength;
    Array<double> j_per;
    Array<double> j_pef;
    Array<double> j_h;

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(ARModeler);
};

#endif AR_MODELER_H_INCLUDED