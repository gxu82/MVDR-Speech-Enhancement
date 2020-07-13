//
// Created by cxk131330 on 10/5/2017.
//

#ifndef AUDIOTWOMICS_2_MVDR2_H
#define AUDIOTWOMICS_2_MVDR2_H


#include <stdio.h>

#include "Transform.h"



typedef struct mvdr_parameters {
    int num_point,num_channel,frame_len,stft_len,frame_shift, framecounter;
    double *sub_data1,*sub_data2,*data1,*data2,*spectrum,*mag,*output,*win,*sys_win,*in_buffer1,*in_buffer2,*in_prev1,*in_prev2,*output_final,*output_old;
    double  *Z_real1,*Z1_real1,*Z_img1,*Z1_img1,*Z_real2,*Z1_real2,*Z_img2,*Z1_img2,*diag1,*diag2;
    int *time,*alpha;

    double energy,sum;

    int frame_count,num_stat,len1, len2,PERC;
    float sum1,determinant;

    double  *z1_real,*z2_real,*z1_img,*z2_img,*rec_signal1_r,*rec_signal2_r,*rec1_r ,*rec2_r,*rec_signal1_i,*rec_signal2_i,*rec1_i,*rec2_i;

    double  *global_covar_r1,*global_covar_r2,*global_covar_r3,*global_covar_r4,*corrvar_r1,*corrvar1_r1,*corrvar_r2,*corrvar1_r2,*corrvar_r3,*corrvar1_r3,*corrvar_r4,*corrvar1_r4,*w_r1,*w_r2,*global_covar_inv_r1,*global_covar_inv_r2,*global_covar_inv_r3,*global_covar_inv_r4,*global_covar_i1,*global_covar_i2,*global_covar_i3,*global_covar_i4,*corrvar_i1,*corrvar1_i1,*corrvar_i2,*corrvar1_i2,*corrvar_i3,*corrvar1_i3,*corrvar_i4,*corrvar1_i4,*w_i1,*w_i2,*global_covar_inv_i1,*global_covar_inv_i2,*global_covar_inv_i3,*global_covar_inv_i4;

    double temp,temp1,temp2,temp3;

    void(*enhance_signal)(struct mvdr_parameters *mvdr, double *input1, double *input2, int framecounter);

} mvdr_parameters;

void enhance_signal(mvdr_parameters *mvdr, double *input1, double *intput2, int framecounter);
mvdr_parameters* newparameters(int stft, int numstat, int framelen, int num_channel, int num_point, int frame_shift);



#endif //AUDIOTWOMICS_2_MVDR2_H
