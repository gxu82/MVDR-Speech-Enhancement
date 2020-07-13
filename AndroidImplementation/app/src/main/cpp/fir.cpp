//
// Created by Shankar, Nikhil on 10/10/17.
//

#include "fir.h"

void fir(double* input, double* output, int nSamples) {
    int i = 0;

    static double endSamples[2] = {0,0};

    for (i = nSamples - 1; i > 1; i--) {
        output[i] = (input[i] + input[i - 1] + input[i - 2])/3;
    }

    output[1] = (input[1] + input[0] + endSamples[0])/3;
    output[0] = (input[0] + endSamples[0] + endSamples[1])/3;

    endSamples[0] = input[nSamples - 1];
    endSamples[1] = input[nSamples - 2];

}
