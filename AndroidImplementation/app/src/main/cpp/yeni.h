//
// Created by axk166230 on 9/7/2017.
//

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdbool.h>

#include <time.h>

#ifndef AUDIOTWOMICS_YENI_H
#define AUDIOTWOMICS_YENI_H

double * filter_m(double * b, double  a, double  *x,  int lenb,  int lena,  int lenx);
double *conv(double *a,double *b,int lena, int lenb);
#endif //AUDIOTWOMICS_YENI_H
