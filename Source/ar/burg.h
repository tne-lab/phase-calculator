// simple header for burg.cpp (AR model estimation using Burg method)

#ifndef BURG_H_INCLUDED
#define BURG_H_INCLUDED

//#define GENERATE_INTERMEDIATES

void ARMaxEntropy (double *inputseries, int length, int degree,
    #ifdef GENERATE_INTERMEDIATES
        double **ar,
    #else
        double *coef,
    #endif
        double *per, double *pef, double *h, double *g);

#endif // BURG_H_INCLUDED