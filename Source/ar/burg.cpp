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
*/

#include "burg.h"

void ARMaxEntropy (double *inputseries, int length, int degree,
    #ifdef GENERATE_INTERMEDIATES
        double **ar,
    #else
        double *coef,
    #endif
        double *per, double *pef, double *h, double *g)
{
    double t1, t2;
    int n;

    for (n = 1; n <= degree; n++)
    {
        double sn = 0.0;
        double sd = 0.0;
        int j;
        int jj = length - n;

        for (j = 0; j < jj; j++)
        {
            t1 = inputseries [j + n] + pef [j];
            t2 = inputseries [j] + per [j];
            sn -= 2.0 * t1 * t2;
            sd += (t1 * t1) + (t2 * t2);
        }

        t1 = g [n] = sn / sd;
        if (n != 1)
        {
            for (j = 1; j < n; j++)
                h [j] = g [j] + t1 * g [n - j];
            for (j = 1; j < n; j++)
                g [j] = h [j];
            jj--;
        }

        for (j = 0; j < jj; j++)
        {
            per [j] += t1 * pef [j] + t1 * inputseries [j + n];
            pef [j] = pef [j + 1] + t1 * per [j + 1] + t1 * inputseries [j + 1];
        }

#ifdef GENERATE_INTERMEDIATES
        for (j = 0; j < n; j++)
           ar [n][j] = g [j + 1];
#endif
    }

#ifndef GENERATE_INTERMEDIATES
    for (n = 0; n < degree; n++)
        coef [n] = g [n + 1];
#endif
}
