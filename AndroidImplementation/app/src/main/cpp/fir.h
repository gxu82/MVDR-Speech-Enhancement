//
// Created by Shankar, Nikhil on 10/10/17.
//

#ifndef AUDIOTWOMICS_2_1_FIR_H
#define AUDIOTWOMICS_2_1_FIR_H

#include "Transform.h"
#include <stdio.h>
void fir(double* input, double* output, int nSamples);
void jmap(double* input, double* output, int nSamples, int Srate, double beta);


#endif //AUDIOTWOMICS_2_1_FIR_H
