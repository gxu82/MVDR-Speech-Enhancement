//
// Created by cxk131330 on 10/5/2017.
//

#include "mvdr2.h"
//#include "mvdr.h"
mvdr_parameters* newparameters(int stft, int numstat, int framelen, int num_channel, int num_point, int frame_shift)
{
    mvdr_parameters* newparameters = (mvdr_parameters*)malloc(sizeof(mvdr_parameters));
    newparameters->framecounter=0;
    newparameters->num_point = num_point;
    newparameters->frame_shift = frame_shift;
    newparameters->num_channel = num_channel;
    newparameters->stft_len = stft;
    newparameters->num_stat = numstat;
    newparameters->frame_len = framelen;

    newparameters->sub_data1 = (double*)calloc(framelen,sizeof(double));
    newparameters->sub_data2 = (double*)calloc(framelen, sizeof(double));
    newparameters->data1 = (double*)calloc(framelen, sizeof(double));
    newparameters->data2 = (double*)calloc(framelen, sizeof(double));
    newparameters->spectrum = (double*)calloc(stft, sizeof(double));
    newparameters->mag = (double*)calloc(stft, sizeof(double));
    newparameters->output = (double*)calloc(framelen,sizeof(double));
    newparameters->Z_real1 = (double*)calloc(stft,sizeof(double));
    newparameters->Z_real2 = (double*)calloc(stft,sizeof(double));
    newparameters->Z_img1 = (double*)calloc(stft,sizeof(double));
    newparameters->Z_img2 = (double*)calloc(stft,sizeof(double));
    newparameters->Z1_real1 = (double*)calloc(stft,sizeof(double));
    newparameters->Z1_real2 = (double*)calloc(stft,sizeof(double));
    newparameters->Z1_img1 = (double*)calloc(stft,sizeof(double));
    newparameters->Z1_img2 = (double*)calloc(stft,sizeof(double));
    newparameters->corrvar_r1 = (double*)calloc(stft/2+1,sizeof(double));
    newparameters->corrvar1_r1 = (double*)calloc(stft/2+1,sizeof(double));
    newparameters->corrvar_r2 = (double*)calloc(stft/2+1,sizeof(double));
    newparameters->corrvar1_r2 = (double*)calloc(stft/2+1,sizeof(double));
    newparameters->corrvar_r3 = (double*)calloc(stft/2+1,sizeof(double));
    newparameters->corrvar1_r3 = (double*)calloc(stft/2+1,sizeof(double));
    newparameters->corrvar_r4 = (double*)calloc(stft/2+1,sizeof(double));
    newparameters->corrvar1_r4 = (double*)calloc(stft/2+1,sizeof(double));

    newparameters->corrvar_i1 = (double*)calloc(stft/2+1,sizeof(double));
    newparameters->corrvar1_i1 = (double*)calloc(stft/2+1,sizeof(double));
    newparameters->corrvar_i2 = (double*)calloc(stft/2+1,sizeof(double));
    newparameters->corrvar1_i2 = (double*)calloc(stft/2+1,sizeof(double));
    newparameters->corrvar_i3 = (double*)calloc(stft/2+1,sizeof(double));
    newparameters->corrvar1_i3 = (double*)calloc(stft/2+1,sizeof(double));
    newparameters->corrvar_i4 = (double*)calloc(stft/2+1,sizeof(double));
    newparameters->corrvar1_i4 = (double*)calloc(stft/2+1,sizeof(double));

    newparameters->global_covar_r1 = (double*)calloc(stft/2+1,sizeof(double));
    newparameters->global_covar_r2 = (double*)calloc(stft/2+1,sizeof(double));
    newparameters->global_covar_r3 = (double*)calloc(stft/2+1,sizeof(double));
    newparameters->global_covar_r4 = (double*)calloc(stft/2+1,sizeof(double));

    newparameters->global_covar_i1= (double*)calloc(stft/2+1,sizeof(double));
    newparameters->global_covar_i2= (double*)calloc(stft/2+1,sizeof(double));
    newparameters->global_covar_i3= (double*)calloc(stft/2+1,sizeof(double));
    newparameters->global_covar_i4= (double*)calloc(stft/2+1,sizeof(double));

    newparameters->w_r1 = (double*)calloc(stft/2+1,sizeof(double));
    newparameters->w_r2 = (double*)calloc(stft/2+1,sizeof(double));

    newparameters->global_covar_inv_r1 = (double*)calloc(stft/2+1,sizeof(double));
    newparameters->global_covar_inv_r2 = (double*)calloc(stft/2+1,sizeof(double));
    newparameters->global_covar_inv_r3 = (double*)calloc(stft/2+1,sizeof(double));
    newparameters->global_covar_inv_r4 = (double*)calloc(stft/2+1,sizeof(double));


    newparameters->w_i1= (double*)calloc(stft/2+1,sizeof(double));
    newparameters->w_i2= (double*)calloc(stft/2+1,sizeof(double));

    newparameters-> global_covar_inv_i1= (double*)calloc(stft/2+1,sizeof(double));
    newparameters-> global_covar_inv_i2= (double*)calloc(stft/2+1,sizeof(double));
    newparameters-> global_covar_inv_i3= (double*)calloc(stft/2+1,sizeof(double));
    newparameters-> global_covar_inv_i4= (double*)calloc(stft/2+1,sizeof(double));

    newparameters->win = (double*)calloc(stft, sizeof(double));
    newparameters->sys_win = (double*)calloc(framelen, sizeof(double));

    newparameters->in_buffer1 = (double*)calloc(framelen,sizeof(double));
    newparameters->in_buffer2 = (double*)calloc(framelen,sizeof(double));
    newparameters->in_prev1 = (double*)calloc(framelen,sizeof(double));
    newparameters->in_prev2 = (double*)calloc(framelen,sizeof(double));

    newparameters->output_final = (double*)calloc(framelen,sizeof(double));
    newparameters->output_old = (double*)calloc(framelen,sizeof(double));

//}
//
    newparameters-> diag1= (double*)calloc(1,sizeof(double));
    newparameters-> diag2= (double*)calloc(1,sizeof(double));
    newparameters->determinant = 0;

    //X3=newTransform(stft_len);
    //X4=newTransform(stft_len);
    newparameters->z1_real= (double*)calloc(stft, sizeof(double));
    newparameters->z2_real= (double*)calloc(stft, sizeof(double));
    newparameters->z1_img= (double*)calloc(stft, sizeof(double));
    newparameters->z2_img= (double*)calloc(stft, sizeof(double));
    //double complex z3[stft_len];
    //double complex z4[stft_len];
    newparameters->rec_signal1_r = (double*)calloc(stft, sizeof(double));
    newparameters->rec_signal2_r= (double*)calloc(stft, sizeof(double));
    newparameters->rec1_r = (double*)calloc(stft, sizeof(double));
    newparameters->rec2_r = (double*)calloc(stft, sizeof(double));

    newparameters->rec_signal1_i = (double*)calloc(stft, sizeof(double));
    newparameters->rec_signal2_i = (double*)calloc(stft, sizeof(double));
    newparameters->rec1_i = (double*)calloc(stft, sizeof(double));
    newparameters->rec2_i = (double*)calloc(stft, sizeof(double));


    newparameters->temp=0;
    newparameters->temp1=0;
    newparameters->temp2=0;
    newparameters->temp3=0;



    newparameters->enhance_signal = enhance_signal;

    return newparameters;

}

void enhance_signal(mvdr_parameters *mvdr, double *input1, double *input2, int framecounter)
{
    Transform *X1;
    Transform *X2;
    Transform *Y;
    //Transform *X3;
    //Transform *X4;
    for (int i = 0; i < mvdr->stft_len; i++)
    {
        mvdr->win[i] = 0.5 * (1 - cosf(2 * M_PI*(i + 1) / (mvdr->stft_len + 1)));
    }
    for (int i = 0; i < mvdr->frame_len; i++)
    {
        mvdr->sys_win[i] = 1 / (mvdr->win[i] + mvdr->win[i + mvdr->frame_len]);
    }
    X1=newTransform(mvdr->stft_len);
    X2=newTransform(mvdr->stft_len);
    Y = newTransform(mvdr->stft_len);
    if (framecounter < mvdr->num_stat)
    {

        for (int i = 0; i < mvdr->frame_len; i++)
        {

           mvdr->sub_data1[i] = input1[i];
            mvdr->sub_data2[i] = input2[i];


            //printf("%d %f\n",i, sub_data1[i]);
        }
        X1->doTransform(X1,mvdr->sub_data1);
        X2->doTransform(X2,mvdr->sub_data2);
        for (int i = 0; i < mvdr->stft_len; i++)
        {
            mvdr->z1_real[i] = X1->real[i];
            mvdr->z1_img[i] = X1->imaginary[i];
            mvdr->z2_real[i] = X2->real[i];
            mvdr->z2_img[i] = X2->imaginary[i];
            mvdr->Z_real1[i]=mvdr->z1_real[i];
            mvdr->Z_img1[i]=mvdr->z1_img[i];
            mvdr->Z_real2[i] = mvdr->z2_real[i];
            mvdr->Z_img2[i] = mvdr->z2_img[i];

        }

        for(int i = 0;i<mvdr->stft_len;i++)
        {

            mvdr->Z1_real1[i] = mvdr->Z_real1[i];
            mvdr->Z1_img1[i] =mvdr-> Z_img1[i]*(-1);
            mvdr->Z1_real2[i] = mvdr->Z_real2[i];
            mvdr->Z1_img2[i] =mvdr-> Z_img2[i]*(-1);
        }


        for(int m = 0; m< (mvdr->stft_len/2)+1;m++)
        {
            mvdr->corrvar_r1[m]=(mvdr->Z_real1[m]*mvdr->Z1_real1[m])-(mvdr->Z_img1[m]*mvdr->Z1_img1[m]);
            mvdr->corrvar_r2[m]=(mvdr->Z_real1[m]*mvdr->Z1_real2[m])-(mvdr->Z_img1[m]*mvdr->Z1_img2[m]);
            mvdr->corrvar_r3[m]=(mvdr->Z_real2[m]*mvdr->Z1_real1[m])-(mvdr->Z_img2[m]*mvdr->Z1_img1[m]);
            mvdr->corrvar_r4[m]=(mvdr->Z_real2[m]*mvdr->Z1_real2[m])-(mvdr->Z_img2[m]*mvdr->Z1_img2[m]);
            mvdr->corrvar1_r1[m]=mvdr->corrvar_r1[m]/(mvdr->corrvar_r1[m]+mvdr->corrvar_r4[m]+8.854*10e-16);
            mvdr->corrvar1_r2[m]=mvdr->corrvar_r2[m]/(mvdr->corrvar_r1[m]+mvdr->corrvar_r4[m]+8.854*10e-16);
            mvdr->corrvar1_r3[m]=mvdr->corrvar_r3[m]/(mvdr->corrvar_r1[m]+mvdr->corrvar_r4[m]+8.854*10e-16);
            mvdr->corrvar1_r4[m] = mvdr->corrvar_r4[m]/(mvdr->corrvar_r1[m]+mvdr->corrvar_r4[m]+8.854*10e-16);

            mvdr->corrvar_i1[m]=(mvdr->Z_real1[m]*mvdr->Z1_img1[m])+(mvdr->Z_img1[m]*mvdr->Z1_real1[m]);
            mvdr->corrvar_i2[m]=(mvdr->Z_real1[m]*mvdr->Z1_img2[m])+(mvdr->Z_img1[m]*mvdr->Z1_real2[m]);
            mvdr->corrvar_i3[m]=(mvdr->Z_real2[m]*mvdr->Z1_img1[m])+(mvdr->Z_img2[m]*mvdr->Z1_real1[m]);
            mvdr->corrvar_i4[m]=(mvdr->Z_real2[m]*mvdr->Z1_img2[m])+(mvdr->Z_img2[m]*mvdr->Z1_real2[m]);
            mvdr->corrvar1_i1[m]=mvdr->corrvar_i1[m]/(mvdr->corrvar_r1[m]+mvdr->corrvar_r4[m]+8.854*10e-16);
            mvdr->corrvar1_i2[m]=mvdr->corrvar_i2[m]/(mvdr->corrvar_r1[m]+mvdr->corrvar_r4[m]+8.854*10e-16);
            mvdr->corrvar1_i3[m]=mvdr->corrvar_i3[m]/(mvdr->corrvar_r1[m]+mvdr->corrvar_r4[m]+8.854*10e-16);
            mvdr->corrvar1_i4[m]=mvdr->corrvar_i4[m]/(mvdr->corrvar_r1[m]+mvdr->corrvar_r4[m]+8.854*10e-16);


        }

//
        for(int m = 0; m< (mvdr->stft_len/2)+1;m++) {
            mvdr->global_covar_r1[m] += mvdr->corrvar1_r1[m];
            mvdr->global_covar_r2[m] += mvdr->corrvar1_r2[m];
            mvdr->global_covar_r3[m] += mvdr->corrvar1_r3[m];
            mvdr->global_covar_r4[m] += mvdr->corrvar1_r4[m];

            mvdr->global_covar_i1[m] += mvdr->corrvar1_i1[m];
            mvdr->global_covar_i2[m] += mvdr->corrvar1_i2[m];
            mvdr->global_covar_i3[m] += mvdr->corrvar1_i3[m];
            mvdr->global_covar_i4[m] += mvdr->corrvar1_i4[m];

        }

    }
    if (framecounter==mvdr->num_stat)
    {
        for (int m = 0; m < (mvdr->stft_len / 2) + 1; m++) {
            mvdr->global_covar_r1[m] = mvdr->global_covar_r1[m] / 129;
            mvdr->global_covar_r2[m] = mvdr->global_covar_r2[m] / 129;
            mvdr->global_covar_r3[m] = mvdr->global_covar_r3[m] / 129;
            mvdr->global_covar_r4[m] = mvdr->global_covar_r4[m] / 129;

            mvdr->global_covar_i1[m] = mvdr->global_covar_i1[m] / 129;
            mvdr->global_covar_i2[m] = mvdr->global_covar_i2[m] / 129;
            mvdr->global_covar_i3[m] = mvdr->global_covar_i3[m] / 129;
            mvdr->global_covar_i4[m] = mvdr->global_covar_i4[m] / 129;
 //__android_log_print(ANDROID_LOG_INFO, "global_covar", "frame= %0.32lf",  mvdr->global_covar_r4[m]);
        }
    }

    if (framecounter >= mvdr->num_stat)
    {

//        if(k+mvdr->frame_len>mvdr->num_point)
//        {
//            break;
//        }
        for (int i = 0; i < mvdr->frame_len; i++)
        {

            mvdr->data1[i] =input1[i];
            mvdr->data2[i] =input2[i];

        }

        X1->doTransform(X1,mvdr->data1);
        X2->doTransform(X2,mvdr->data2);

        for (int i = 0; i < mvdr->stft_len; i++)
        {
            mvdr->z1_real[i] = X1->real[i];
            mvdr->z1_img[i] = X1->imaginary[i];
            mvdr->z2_real[i] = X2->real[i];
            mvdr->z2_img[i] = X2->imaginary[i];
            //printf("%f  i%f\n", creal(z4[i]), cimag(z4[i]));
        }

        for(int k1 =0;k1<mvdr->stft_len/2+1;k1++)
        {
            mvdr->global_covar_r1[k1] = mvdr->global_covar_r1[k1] + 1*1e-08;
            mvdr->global_covar_r4[k1] = mvdr->global_covar_r4[k1] ;

            mvdr->determinant = (mvdr->global_covar_r1[k1]*mvdr->global_covar_r4[k1])-((mvdr->global_covar_r2[k1]*mvdr->global_covar_r3[k1])-(mvdr->global_covar_i2[k1]*mvdr->global_covar_i3[k1]));


            mvdr->global_covar_inv_r1[k1]=mvdr->global_covar_r4[k1]/mvdr->determinant;
            mvdr->global_covar_inv_r2[k1]=-mvdr->global_covar_r2[k1]/mvdr->determinant;
            mvdr->global_covar_inv_r3[k1]=-mvdr->global_covar_r3[k1]/mvdr->determinant;
            mvdr->global_covar_inv_r4[k1]=mvdr->global_covar_r1[k1]/mvdr->determinant;

            mvdr->global_covar_inv_i1[k1]=mvdr->global_covar_i4[k1]/mvdr->determinant;
            mvdr->global_covar_inv_i2[k1]=-mvdr->global_covar_i2[k1]/mvdr->determinant;
            mvdr->global_covar_inv_i3[k1]=-mvdr->global_covar_i3[k1]/mvdr->determinant;
            mvdr->global_covar_inv_i4[k1]=mvdr->global_covar_i1[k1]/mvdr->determinant;




            mvdr->w_r1[k1] = (mvdr->global_covar_inv_r1[k1] + mvdr->global_covar_inv_r2[k1])/ (mvdr->global_covar_inv_r1[k1] + mvdr->global_covar_inv_r2[k1] + mvdr->global_covar_inv_r3[k1] + mvdr->global_covar_inv_r4[k1]);
            mvdr->w_r2[k1] = (mvdr->global_covar_inv_r3[k1] + mvdr->global_covar_inv_r4[k1])/(mvdr->global_covar_inv_r1[k1] + mvdr->global_covar_inv_r2[k1] + mvdr->global_covar_inv_r3[k1] + mvdr->global_covar_inv_r4[k1]);

            mvdr->w_i1[k1] = (mvdr->global_covar_inv_i1[k1] + mvdr->global_covar_inv_i2[k1])/ (mvdr->global_covar_inv_r1[k1] + mvdr->global_covar_inv_r2[k1] + mvdr->global_covar_inv_r3[k1] + mvdr->global_covar_inv_r4[k1]);
            mvdr->w_i2[k1] = (mvdr->global_covar_inv_i3[k1] + mvdr->global_covar_inv_i4[k1])/(mvdr->global_covar_inv_r1[k1] + mvdr->global_covar_inv_r2[k1] + mvdr->global_covar_inv_r3[k1] + mvdr->global_covar_inv_r4[k1]);


        }
        for(int i = 0; i < mvdr->stft_len/2 +1; i++)
        {
            mvdr->w_i1[i] = mvdr->w_i1[i]*(-1);
            mvdr->w_i2[i] = mvdr->w_i2[i]*(-1);

            mvdr->rec_signal1_r[i] = (mvdr->w_r1[i]*mvdr->z1_real[i])-(mvdr->w_i1[i]*mvdr->z1_img[i]);
            mvdr->rec_signal1_i[i] = (mvdr->w_r1[i]*mvdr->z1_img[i])+(mvdr->w_i1[i]*mvdr->z1_real[i]);
            // rec_signal2[i] = w[1][i]*z2[i];
            mvdr->rec_signal2_r[i] = (mvdr->w_r2[i]*mvdr->z2_real[i])-(mvdr->w_i2[i]*mvdr->z2_img[i]);
            mvdr->rec_signal2_i[i] = (mvdr->w_r2[i]*mvdr->z2_img[i])+(mvdr->w_i2[i]*mvdr->z2_real[i]);


        }
        for(int i = 1; i < mvdr->stft_len/2; i++)
        {
            mvdr->rec1_r[i-1] = mvdr->rec_signal1_r[i];
            mvdr->rec2_r[i-1] = mvdr->rec_signal2_r[i];

            mvdr->rec1_i[i-1] = mvdr->rec_signal1_i[i];
            mvdr->rec2_i[i-1] = mvdr->rec_signal2_i[i];
        }

        int i =mvdr->stft_len/2 -1;
        int j = i - 1;   // j will Point to last Element
        i = 0;       // i will be pointing to first element

        while (i < j) {
            mvdr->temp = mvdr->rec1_r[i];
            mvdr->rec1_r[i] = mvdr->rec1_r[j];
            mvdr->rec1_r[j] = mvdr->temp;
            i++;             // increment i
            j--;          // decrement j
        }
        i =mvdr->stft_len/2 -1;
        j = i - 1;   // j will Point to last Element
        i = 0;       // i will be pointing to first element
        while (i < j) {
            mvdr->temp1 = mvdr->rec2_r[i];
            mvdr->rec2_r[i] = mvdr->rec2_r[j];
            mvdr->rec2_r[j] = mvdr->temp1;
            i++;             // increment i
            j--;          // decrement j
        }
        i =mvdr->stft_len/2 -1;
        j = i - 1;   // j will Point to last Element
        i = 0;       // i will be pointing to first element
        while (i < j) {
            mvdr->temp2 = mvdr->rec1_i[i];
            mvdr->rec1_i[i] = mvdr->rec1_i[j];
            mvdr->rec1_i[j] = mvdr->temp2;
            i++;             // increment i
            j--;          // decrement j
        }

        i =mvdr->stft_len/2 -1;
        j = i - 1;   // j will Point to last Element
        i = 0;       // i will be pointing to first element
        while (i < j) {
            mvdr->temp3 = mvdr->rec2_i[i];
            mvdr->rec2_i[i] = mvdr->rec2_i[j];
            mvdr->rec2_i[j] = mvdr->temp3;
            i++;             // increment i
            j--;          // decrement j
        }
        for(int i = 0; i < mvdr->stft_len/2-1 ; i++)
        {
            mvdr->rec1_i[i] = (-1)*(mvdr->rec1_i[i]);
            mvdr->rec2_i[i] = (-1)*(mvdr->rec2_i[i]);
            // printf("The product: Z1 x Z2 = %f + i%f %d\n", creal(rec_signal1[i]), cimag(rec_signal1[i]), i);
        }

        for(int i = 0; i< mvdr->stft_len/2-1; i++)
        {
            //rec_signal1[i] = rec_signal[i];
            mvdr->rec_signal1_r[i+mvdr->stft_len/2+1] = mvdr->rec1_r[i];
            mvdr->rec_signal2_r[i+mvdr->stft_len/2+1] = mvdr->rec2_r[i];

            mvdr->rec_signal1_i[i+mvdr->stft_len/2+1] = mvdr->rec1_i[i];
            mvdr->rec_signal2_i[i+mvdr->stft_len/2+1] = mvdr->rec2_i[i];
        }

        for(int i =0; i<mvdr->stft_len;i++)
        {
            mvdr->rec_signal1_r[i] = mvdr->rec_signal1_r[i] + mvdr->rec_signal2_r[i];
            mvdr->rec_signal1_i[i] = mvdr->rec_signal1_i[i] + mvdr->rec_signal2_i[i];
            //printf("%d %f\n",i, (rec_signal1_i[i]));
            X1->real[i] = mvdr->rec_signal1_r[i];
            X1->imaginary[i] = mvdr->rec_signal1_i[i];

            // printf("The product: Z1 x Z2 = %f + i%f %d\n", creal(rec_signal1[i]), cimag(rec_signal1[i]), i);
        }
        Y->invTransform(Y,X1->real,X1->imaginary);

        for(int i =0; i<mvdr->frame_len;i++)
        {
            mvdr->output[i] = Y->real[i];
//            mvdr->output_final[i]=(mvdr->output_old[i]+Y->real[i])*mvdr->sys_win[i];
//            mvdr->output_old[i]=Y->real[i+mvdr->frame_len];

        }
        for(int i =0; i<mvdr->frame_len;i++) {
            //mvdr->output[i] = mvdr->output_final[i];
        }
    }
    if(framecounter<mvdr->num_stat)
    {
        for (int i = 0; i < mvdr->frame_len; i++) {


            mvdr->output[i] = input1[i];
        }
    }

}

