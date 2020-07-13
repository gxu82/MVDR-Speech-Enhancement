//
// Created by Abdullah on 9/4/2017.
//

#ifndef AUDIOTWOMICS_TRANSFORM_H
#define AUDIOTWOMICS_TRANSFORM_H

#define _USE_MATH_DEFINES
#define PI 3.14159265
#define MAXIT 1000
#define EULER 0.5772156649
#define FPMIN 1.0e-30
#define EPSS  1.0e-7f
#define pi 3.1415926535897932384626433832795

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include <android/log.h>
#include <time.h>

typedef struct Transform {
    int points;
    double* sine;
    double* cosine;
    double* real;
    double* imaginary;
    void(*doTransform)(struct Transform* transform, double* input);
    void(*invTransform)(struct Transform* transform, double* inputreal, double* inputimaginary);
} Transform;

Transform* newTransform(int points);
void transformMagnitude(Transform* transform, double* output);
void invtranMagnitude(Transform* transform, double* output);
void destroyTransform(Transform* transform);
void FFT(Transform* fft, double* input);
void IFFT(Transform* fft, double* inputreal, double* inputimaginary);

#endif //AUDIOTWOMICS_TRANSFORM_H
