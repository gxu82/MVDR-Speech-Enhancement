//
// Created by Abdullah Kucuk on 9/7/2017.
//

#include "yeni.h"




double * filter_m(double * b, double  a, double  *x,  int lenb,  int lena,  int lenx)// modified filter, it filters when a=1
{

    /*for (int i = 0; i < 223; i++)
    {
    //printf("%lf ", norm_freq[i]);

    //printf("%fi ", x);
    //printf("emphasis \n",i );
    //for (int k=0;k<3;k++)
    printf("fir2  %1.7g \t i=%d	\n", b[i], i);
    }*/
    double * y;
    unsigned int n, i;
    double temp1, temp2;
    y = (double *)malloc(lenx * sizeof(double));

    for (n = 0; n<lenx; n++)
    {
        temp1 = 0;
        temp2 = 0;

        //FeedForward Section Computation
        for (i = 0; (i <= n) && (i<lenb); i++)
            temp1 = temp1 + b[i] * x[(n - i)];

        //Recursive Section Computation
        for (i = 1; (i <= n) && (i<lena); i++)
            temp2 = temp2 + y[n - i];

        y[n] = (temp1 - temp2);


    }
    return y;
}

double *conv(double *a,double *b,int lena, int lenb)
{
    double *c =(double*) calloc((lena+lenb-1),sizeof(double));
    int counter = 1;
    while (counter<lena+lenb)
    {
        if(counter<lenb)
        {
            for(int k=0;k<counter;k++)
            {
                c[counter-1] += a[counter-k-1]*b[k];
            }
        }
        else if(counter<=lena)
        {
            for(int k=0;k<lenb;k++)
            {
                c[counter-1] += a[counter-k-1]*b[k];
            }
        }
        else //tail part
        {
            for(int k=0;k<lena+lenb-counter;k++)
            {
                c[counter-1] += a[lena-k-1]*b[counter-lena+k];
            }
        }
        counter++;
    }
    return c;
}