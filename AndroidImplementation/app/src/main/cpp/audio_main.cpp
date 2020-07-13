//
// Created by Abdullah Kucuk on 9/4/2017.
//


#include <cassert>
#include <cstring>
#include <jni.h>

#include <sys/types.h>
#include <SLES/OpenSLES.h>

#include "audio_common.h"
//#include "audio_recorder.h"
//#include "audio_player.h"
#include "audio_main.h"
#include "mvdr2.h"
//#include "twoDOA.h"

struct EchoAudioEngine {
    SLmilliHertz fastPathSampleRate_;
    uint32_t     fastPathFramesPerBuf_;
    uint16_t     sampleChannels_;
    uint16_t     bitsPerSample_;

    SLObjectItf  slEngineObj_;
    SLEngineItf  slEngineItf_;

    /* AudioRecorder  *recorder_;
     AudioPlayer    *player_;
     AudioQueue     *freeBufQueue_;    //Owner of the queue
     AudioQueue     *recBufQueue_;     //Owner of the queue

     sample_buf  *bufs_;*/
    uint32_t     bufCount_;
    uint32_t     frameCount_;
};
static EchoAudioEngine engine;


bool EngineService(void* ctx, uint32_t msg, void* data );

extern "C" {
//JNIEXPORT void JNICALL
//Java_com_google_sample_echo_MainActivity_createSLEngine(JNIEnv *env, jclass, jint, jint);
JNIEXPORT jlong JNICALL
        Java_com_example_abdullah_audiotwomics_SignalProcessing_paramInitilization(JNIEnv *env,
                                                                                   jobject thiz,
                                                                                   int frequency,
                                                                                   int stepSize,
                                                                                   int frameSize,
                                                                                   float stepFactor,
                                                                                   int noisetype,
                                                                                   bool isEnchanced,
                                                                                   bool isOn);
JNIEXPORT void JNICALL
        Java_com_example_abdullah_audiotwomics_SignalProcessing_realtimeProcessing(JNIEnv *env,
                                                                                   jobject thiz,
                                                                                   jlong memoryPointer,
                                                                                   jshortArray input);
JNIEXPORT void JNICALL
        Java_com_example_abdullah_audiotwomics_SignalProcessing_paramElimination(JNIEnv *env, jobject thiz,
                                                                                 jlong memoryPointer);
JNIEXPORT jshortArray JNICALL
        Java_com_example_abdullah_audiotwomics_SignalProcessing_soundOutput(JNIEnv *env, jobject thiz,
                                                                            jlong memoryPointer,
                                                                            jint outputSelection);

JNIEXPORT jfloatArray JNICALL
        Java_com_example_abdullah_audiotwomics_SignalProcessing_FFTtiming(JNIEnv *env, jobject thiz,
                                                                          jlong memoryPointer,
                                                                          jint outputSelection);

JNIEXPORT jshortArray JNICALL
        Java_com_example_abdullah_audiotwomics_SignalProcessing_dataOutput(JNIEnv *env, jobject thiz,
                                                                           jlong memoryPointer,
                                                                           jint outputSelection);
JNIEXPORT jfloat JNICALL
        Java_com_example_abdullah_audiotwomics_SignalProcessing_getTime(JNIEnv *env, jobject thiz,
                                                                        jlong memoryPointer);
JNIEXPORT jfloat JNICALL
        Java_com_example_abdullah_audiotwomics_SignalProcessing_getComputeTime(JNIEnv *env, jobject thiz,
                                                                               jlong memoryPointer);
JNIEXPORT jfloat JNICALL
        Java_com_example_abdullah_audiotwomics_SignalProcessing_getFilteringTime(JNIEnv *env, jobject thiz,
                                                                                 jlong memoryPointer);





/*
JNIEXPORT jboolean JNICALL
Java_com_google_sample_echo_MainActivity_createSLBufferQueueAudioPlayer(JNIEnv *env, jclass);
JNIEXPORT void JNICALL
Java_com_google_sample_echo_MainActivity_deleteSLBufferQueueAudioPlayer(JNIEnv *env, jclass type);

JNIEXPORT void JNICALL
Java_com_google_sample_echo_SignalProcessing_realtimeProcessing(JNIEnv *env, jobject thiz,  jlong memoryPointer, jshortArray input);
JNIEXPORT void JNICALL
Java_com_google_sample_echo_MainActivity_deleteAudioRecorder(JNIEnv *env, jclass type);
JNIEXPORT void JNICALL

}*/
JNIEXPORT void JNICALL
        Java_com_example_abdullah_audiotwomics_SignalProcessing_realtimeProcessing(JNIEnv *env,
                                                                                   jobject thiz,
                                                                                   jlong memoryPointer,
                                                                                   jshortArray input);

JNIEXPORT void JNICALL
        Java_com_example_abdullah_audiotwomics_MainActivity_startPlay(JNIEnv *env, jclass type);
JNIEXPORT void JNICALL
        Java_com_example_abdullah_audiotwomics_MainActivity_stopPlay(JNIEnv *env, jclass type);


}

JNIEXPORT jlong JNICALL
Java_com_example_abdullah_audiotwomics_SignalProcessing_paramInitilization(
        JNIEnv *env, jobject thiz,int frequency, int stepSize,int frameSize, float stepFactor, int noisetype, bool isEndFire, bool isOn){
    Variables *mainParam = (Variables*)malloc(sizeof(Variables));
    int num_point = 160000;
    int num_channel = 2;
    int frame_len = 960;
    int stft_len = 2048;
    int frame_shift = 480;
    int num_stat = 129;

    int i;
    int beam_filtlen=2048;
    mainParam->counter = 0;
    mainParam->u = stepFactor;
    mainParam->stepSize = stepSize;
    mainParam->samplingRate = frequency;
    mainParam->frameSize = frameSize;  ////Threadshold Time
    mainParam->filLen = 1;
    mainParam->topMicBuffer = (double*)calloc(frameSize+mainParam->filLen-1,sizeof(double));
    mainParam->botMicBuffer = (double*)calloc(frameSize+mainParam->filLen-1,sizeof(double));
    mainParam->originalBuffer = (int*)calloc(2*stepSize,sizeof(int));
    mainParam->mixedBuffer = (double*)calloc(stepSize,sizeof(double));
    mainParam->w = (float*)malloc(mainParam->filLen*sizeof(float));
    mainParam->filter_output1 = (double*)calloc(stepSize+beam_filtlen,sizeof(double));
    mainParam->filter_output2 = (double*)calloc(stepSize+beam_filtlen,sizeof(double));
    mainParam->output = (double*)calloc(beam_filtlen,sizeof(double));
    mainParam->prev_filter_output1 = (double*)calloc(stepSize,sizeof(double));
    mainParam->prev_filter_output2 = (double*)calloc(stepSize,sizeof(double));

    for(i=0;i<mainParam->filLen;i++){
        mainParam->w[i]=1;
    }
    mainParam->y_curr = (float*)calloc(frameSize,sizeof(float));
    mainParam->y_prev = (float*)calloc(frameSize,sizeof(float));
    mainParam->y = (float*)calloc(stepSize,sizeof(float));
    mainParam->e = (float*)calloc(stepSize,sizeof(float));
    mainParam->outputBuffer = (double*)calloc(stepSize,sizeof(double));
    mainParam->uplimit = 5*frequency/stepSize;
    //mainParam->logmmsePtr = initiallogMMSE(frameSize, stepSize, frequency,noisetype);
    //mainParam->timer = newTimer();
    /*mainParam->Threadtime = ThreadTime;
    mainParam->DurationTime = DurationTime;
    mainParam->angle = 0;*/
    mainParam->angle_count = 0;
    // mainParam->trans1 = newTransform(DOA_NFFT);
    //mainParam->trans2 = newTransform(DOA_NFFT);
    mainParam->isEndFire=isEndFire;
    mainParam->isOn=isOn;
    //mainParam->logMMSE=newparameters(DOA_NFFT, DOA_n, 50); //PERC=50
    //mainParam->logMMSE2=newparameters(DOA_NFFT, DOA_n, 50); //PERC=50
    mainParam->mvdr = newparameters(stft_len,num_stat,frame_len,num_channel,num_point,frame_shift);
    mainParam->ake_counter=0;
    mainParam->ake_avg_timer=0;

    __android_log_print(ANDROID_LOG_ERROR, "Parameter Initialized","Successfully");
    return (jlong)mainParam;
}


JNIEXPORT void JNICALL
Java_com_example_abdullah_audiotwomics_SignalProcessing_realtimeProcessing(JNIEnv *env, jobject thiz,  jlong memoryPointer, jshortArray input){
    Variables *mainParam = (Variables*) memoryPointer;
    //startTimer(mainParam->timer);
    short *_in = env->GetShortArrayElements(input, NULL);
    int i,j,stepSize,frameSize,filLen;
    float tempVar,tempEng, mThreadTime, mDurationTime;
    stepSize = mainParam->stepSize;
    frameSize = mainParam->frameSize;
    filLen = mainParam->filLen;
    //mThreadTime = mainParam->Threadtime;
    //mDurationTime = mainParam->DurationTime;
    int beam_filtlen=512;

    double * x_frame = (double*)malloc(mainParam->stepSize * sizeof(double));
    double * h_frame = (double*)malloc(mainParam->stepSize * sizeof(double));


    for (int k = 0; k < mainParam->stepSize; k++) {

        x_frame[k] = mainParam->botMicBuffer[k];
        h_frame[k] = mainParam->topMicBuffer[k];
        //  __android_log_print(ANDROID_LOG_INFO, "global_covar", "h_frame= %1.9lf", (h_frame[k]));

    }


    clock_t t;
    t = clock();

    mainParam->ake_counter= mainParam->ake_counter+1;

    if(mainParam->isOn) {

        for (int i = 0; i < nFFT; i++)
        {
            mainParam->mvdr->win[i] = 0.5 * (1 - cosf(2 * M_PI*(i + 1) / (nFFT + 1)));
        }
        for (int i = 0; i < FRAMESIZE; i++)
        {
            mainParam->mvdr->sys_win[i] = 1 / (mainParam->mvdr->win[i] + mainParam->mvdr->win[i + FRAMESIZE]);
        }






        mainParam->mvdr->enhance_signal(mainParam->mvdr, x_frame, h_frame, mainParam->ake_counter);

        for (int i = 0; i < FRAMESIZE; i++) {
           //input1[i] = mainParam->mvdr->output[i];
            input1[i] = mainParam->mvdr->output[i];
          //  __android_log_print(ANDROID_LOG_INFO, "global_covar", "frame= %.32lf",   input1[i]);
            // printf("%f\n",input[i]);
        }


        X=newTransform(nFFT);
        Y=newTransform(nFFT);

        X->doTransform(X,in_buffer);
        transformMagnitude(X, fftmagoutput);


       // mainParam->ake_counter++;
        for(int i =0;i<FRAMESIZE;i++)
        {
            in_buffer[i]=in_prev[i] * mainParam->mvdr->win[i];
            in_prev[i]=input1[i];
            in_buffer[i+FRAMESIZE]=input1[i] * mainParam->mvdr->win[i+FRAMESIZE];
            //__android_log_print(ANDROID_LOG_INFO, "global_covar", "frame= %.32lf",   mainParam->mvdr->win[i]);

        }
        float total_noisepower=0;
        float total_speechpower=0;
        float sum_ensig = 0;
        float sum_nsig = 0;


        if(mainParam->ake_counter<=6)
        {
            for (int i = 0; i < nFFT; i++)
            {
                noise_mean[i] += fftmagoutput[i];

                //printf("%f\n",noise_mean[i]);
            }
        }

        if(mainParam->ake_counter==7)
        {
            for (int i = 0; i < nFFT; i++)
            {
                noise_mu2[i]=pow(noise_mean[i]/7,2);

                //printf("%f\n",noise_mu2[i]);

            }
        }
        if(mainParam->ake_counter>=7)
        {
            for (int i = 0; i < nFFT; i++)
            {
                sig2[i] = pow(fftmagoutput[i], 2);


                if( noise_mu2[i]==0)
                {
                    noise_mu2[i] = noise_mu2[i] + epsilon;
                }
                gammak[i] = sig2[i] / noise_mu2[i] < 40 ? sig2[i] / noise_mu2[i] : 40;
                // printf("%f\n",gammak[i]);

                if (mainParam->ake_counter == 7)
                {
                    // for (int i = 0; i < nFFT; i++) {
                    max = gammak[i] - 1 > 0 ? gammak[i] - 1 : 0;
                    ksi[i] = aa + (1 - aa) * max;
                    // printf("%f\n",ksi[i]);
                    // }
                }
                else
                {

                    max = gammak[i] - 1 > 0 ? gammak[i] - 1 : 0;
                    ksi[i] = aa * (Xk_prev[i] / noise_mu2[i]) + (1 - aa) * max;
                    ksi[i] = ksi_min > ksi[i] ? ksi_min : ksi[i];
                    // printf("%f\n",Xk_prev[i]);
                }


                sum_log_sigma_k = 0;
                log_sigma_k[i] = gammak[i] * ksi[i] / (1 + ksi[i]) - log(1 + ksi[i]);
                // if(log_sigma_k[i]<100)
                // {
                sum_log_sigma_k += log_sigma_k[i];
                //
            }
            // printf("%f\n",sum_log_sigma_k);
            vad_decision = sum_log_sigma_k / (nFFT);

            //__android_log_print(ANDROID_LOG_VERBOSE, "MyApp", "The value of SNR_db is %f",vad_decision );
            if (vad_decision < eta)// % noise on
            {
                count = count + 1;
                //printf("%f\n",ksi[i]);
                for (int i = 0; i < nFFT; i++)
                {
                    noise_pow[i] = noise_pow[i] + noise_mu2[i];
                    noise_mu2[i] = noise_pow[i] / count;
                }
            }
            for (int i = 0; i < nFFT; i++)
            {
                vk[i] = ksi[i] * gammak[i] / (1 + ksi[i]);

                hw[i] = (ksi[i] + sqrt(pow(ksi[i], 2) + ((2*beta)*(beta + ksi[i])) * ksi[i] / (gammak[i]+epsilon))) / (2*beta * (beta + ksi[i]));

                sum_ensig += pow(ensig[i], 2);
   
                sum_nsig += pow(fftmagoutput[i], 2);
           
            }

            float PR = sum_ensig/sum_nsig;
            // printf("%f\n",PR);

            if (PR >= 0.4)
                PRT = 1;
            else
                PRT = PR;

            if (PRT == 1)
                N = 1;
            else
                N = 2 * round((1 - PRT / 0.4) * 10 ) + 1;

            //printf("%f\n",N);
            for (int i = 0; i < (int) N; i++)
                H[i] = 1 / N;
            //	H(1:N) = 1 / N;

            for (int i = 0; i < N + nFFT - 1; i++)
            {
                int kmin, kmax, k;

                HPF[i] = 0;

                kmin = (i >= nFFT - 1) ? i - (nFFT - 1) : 0;
                kmax = (i < N - 1) ? i : N - 1;

                for (k = kmin; k <= kmax; k++)
                    HPF[i] += H[k] * sqrt(hw[i - k]* hw[i - k]);
            }
            for (int i = 0; i < nFFT; i++)
            {
                fftmagoutput[i] = fftmagoutput[i] * HPF[i];

                X->real[i]=X->real[i]*HPF[i];
                X->imaginary[i]=X->imaginary[i]*HPF[i];
                // printf("%f\n",X->real[i]);
                //  printf("%f\n",X->imaginary[i]);

                //fftmagoutput[i]=sqrt( X->real[i] * X->real[i] + X->imaginary[i] * X->imaginary[i]);

                ensig[i]=fftmagoutput[i];
                //printf("%f\n",ensig[i]);

                //self->magnitudeRight[i] = fftmagoutput[i];
                if (fftmagoutput[i]<100)
                {
                    // total_speechpower+=magnitudeLeft[i]*magnitudeLeft[i];
                    total_speechpower+=pow(fftmagoutput[i],2);
                }
                if(noise_pow[i]<100)
                {
                    total_noisepower+=noise_pow[i];
                }

                //  }
                // for (int i = 0; i < nFFT; i++)
                // {
                Xk_prev[i] = pow(fftmagoutput[i],2);
            }
        }
        sum_SNR+=SNR_db;
        snr_counter++;
        if (snr_counter==100)
        {
            //pthread_mutex_lock(&mutex);
            snr_avg=sum_SNR/snr_counter;
            //  __android_log_print(ANDROID_LOG_VERBOSE, "MyApp", "The value of SNR_db is %f",snr_counter );
            sum_SNR=0;
            snr_counter=0;
        }
        Y->invTransform(Y,X->real,X->imaginary);


        for(int i=0;i<FRAMESIZE;i++)
        {
            output_final[i]=(output_old[i]+Y->real[i])*mainParam->mvdr->sys_win[i];
            output_old[i]=Y->real[i+FRAMESIZE];
        }




        fir(output_final,float_output,FRAMESIZE);

            //fir(Y->real,float_output,FRAMESIZE);




            for(int i =0;i<mainParam->mvdr->frame_len;i++)
            {
                (mainParam->output[i]) = output_final[i];


            }

        }


    t = clock() - t;

    float ake_timing=(((float)t)/ CLOCKS_PER_SEC)+(mainParam->ake_avg_timer*(mainParam->ake_counter));


    mainParam->ake_avg_timer=ake_timing/(mainParam->ake_counter);

    mainParam->angle_count++;

    for(i = 0;i<filLen-1;i++){
        mainParam->topMicBuffer[i] = mainParam->topMicBuffer[i+stepSize];
        mainParam->botMicBuffer[i] = mainParam->botMicBuffer[i+stepSize];
    }
    for(i = filLen-1,j=0;i<filLen+stepSize-1;i++,j+=2){
        mainParam->topMicBuffer[i]= mainParam->topMicBuffer[i+stepSize];
        mainParam->topMicBuffer[i+stepSize]= _in[j+1]*FAC;
        //mainParam->topMicBuffer[i+stepSize]= _in[2*j]*FAC;
        //mainParam->topMicBuffer[i+stepSize+1]= _in[2*j+1]*FAC;
        mainParam->botMicBuffer[i]= mainParam->botMicBuffer[i+stepSize];
        mainParam->botMicBuffer[i+stepSize]= _in[j]*FAC;
        // mainParam->botMicBuffer[i+stepSize]= _in[2*j+2]*FAC;
        //mainParam->botMicBuffer[i+stepSize+1]= _in[2*j+3]*FAC;
        mainParam->originalBuffer[j] = _in[j];
        mainParam->originalBuffer[j+1] = _in[j+1];
    }

    if(mainParam->isOn){
        for(i=0;i<stepSize;i++){

            mainParam->mixedBuffer[i] = (mainParam->output[i]);
        }
    }else{
        for(i=0;i<stepSize;i++){
            mainParam->mixedBuffer[i] = (x_frame[i]+h_frame[i])/2;
        }

    }

    //  __android_log_print(ANDROID_LOG_ERROR, "Parameter Computed 1st","Successfully");
    env->ReleaseShortArrayElements(input, _in, 0);
    tempEng = 0;

    // *************************************************************************************************** //

    // *************************************************************************************************** //
    // startfilteringTimer(mainParam->timer);
    for(i=0;i<frameSize;i++){
        mainParam->y_prev[i] = mainParam->y_curr[i];
        tempVar = 0;
        for(j=0;j<filLen;j++){
            tempVar += mainParam->w[j]*mainParam->botMicBuffer[i+j];
            tempEng += mainParam->botMicBuffer[i+j]*mainParam->botMicBuffer[i+j];
        }
        mainParam->y_curr[i] = tempVar;
    }


    for(i=0;i<stepSize;i++){
        mainParam->e[i] = mainParam->topMicBuffer[i+filLen-1]-0.5*(mainParam->y_prev[i+stepSize]+mainParam->y_curr[i]);
    }
    //__android_log_print(ANDROID_LOG_ERROR, "Parameter Computed 2nd","Successfully");
    for(i=0;i<filLen;i++){
        tempVar = 0;
        for(j=0;j<stepSize;j++){
            tempVar += mainParam->botMicBuffer[filLen-1+j]*mainParam->e[j];
            tempVar += mainParam->botMicBuffer[j+i]*mainParam->e[j];
        }
        mainParam->w[i] += mainParam->u*tempVar/(tempEng+EPS);
    }
    // stopfilteringTimer(mainParam->timer);
    // *************************************************************************************************** //

    // *************************************************************************************************** //

    if(mainParam->counter<mainParam->uplimit){
        mainParam->counter += 1;
        for(i = 0;i<stepSize;i++){
            mainParam->outputBuffer[i] = mainParam->e[i];
        }
    }else{
        // startcomputeTimer(mainParam->timer);
        //logMMSEexe(mainParam->logmmsePtr, mainParam->e, mainParam->outputBuffer);
        // stopcomputeTimer(mainParam->timer);
    }

    // *************************************************************************************************** //

    // *************************************************************************************************** //
    // stopTimer(mainParam->timer);
}


JNIEXPORT void JNICALL
Java_com_example_abdullah_audiotwomics_SignalProcessing_paramElimination(JNIEnv* env, jobject thiz, jlong memoryPointer){
    Variables *mainParam = (Variables*) memoryPointer;
    if(mainParam != NULL){
        free(mainParam->topMicBuffer);mainParam->topMicBuffer = NULL;
        free(mainParam->botMicBuffer);mainParam->botMicBuffer = NULL;
        free(mainParam->outputBuffer);mainParam->outputBuffer = NULL;
        free(mainParam->originalBuffer);mainParam->originalBuffer = NULL;
        free(mainParam->mixedBuffer);mainParam->mixedBuffer = NULL;
        free(mainParam->y_prev);mainParam->y_prev = NULL;
        free(mainParam->y_curr);mainParam->y_curr = NULL;
        free(mainParam->y);mainParam->y = NULL;
        free(mainParam->e);mainParam->e = NULL;
        free(mainParam->w);mainParam->w = NULL;
        //destroylogMMSE(&(mainParam->logmmsePtr));
        // destroyTimer(&(mainParam->timer));
        free(mainParam);mainParam = NULL;
    }
}

JNIEXPORT jshortArray JNICALL
Java_com_example_abdullah_audiotwomics_SignalProcessing_soundOutput(JNIEnv* env, jobject thiz, jlong memoryPointer, jint outputSelection){
    Variables* mainParam = (Variables*) memoryPointer;
    jshortArray output = env->NewShortArray(mainParam->stepSize);
    short *_output = env->GetShortArrayElements( output, NULL);
    int i;
    switch (outputSelection){
        case 0:		// Original

            for (i=0;i<mainParam->stepSize;i++){
                _output[i] = (short)checkRange(mainParam->mixedBuffer[i]);
            }

            break;
        case 1:		// Enhanced

            break;
        case 2:
            break;
    }
    env->ReleaseShortArrayElements(output, _output, 0);
    return output;
}


JNIEXPORT jfloatArray JNICALL
Java_com_example_abdullah_audiotwomics_SignalProcessing_FFTtiming(JNIEnv* env, jobject thiz, jlong memoryPointer, jint outputSelection){
    Variables* mainParam = (Variables*) memoryPointer;
    jfloatArray output = env->NewFloatArray(2);
    //float tmp = 0;
    float *_output = env->GetFloatArrayElements( output, NULL);



    //int i;
    switch (outputSelection){
        case 0:		// Original



            // tmp = max1(2,3);

            // mainParam->angle = tmp;
            // mainParam->angle_count++;
            _output[0] = mainParam->ake_counter;
            _output[1] = mainParam->ake_avg_timer;

            break;
        case 1:		// Enhanced
            _output[0] = -1;

            break;
        case 2:
            break;
    }
    env->ReleaseFloatArrayElements(output, _output, 0);
    return output;
}

JNIEXPORT jshortArray JNICALL
Java_com_example_abdullah_audiotwomics_SignalProcessing_dataOutput(JNIEnv* env, jobject thiz, jlong memoryPointer, jint outputSelection){
    Variables* mainParam = (Variables*) memoryPointer;
    jshortArray output = env->NewShortArray(2*mainParam->stepSize);
    short *_output = env->GetShortArrayElements(output, NULL);
    int i;
    switch (outputSelection){
        case 0:
            for(i= 0;i<2*mainParam->stepSize;i++){
                _output[i] = mainParam->originalBuffer[i];
            }
            break;
        case 1:
            for(i= 0;i<2*mainParam->stepSize;i++){
                _output[i] = (short)checkRange(mainParam->outputBuffer[i]);
            }
            break;
        case 2:
            break;
    }
    env->ReleaseShortArrayElements(output, _output, 0);
    return output;
}






short checkRange(double input)
{
    int output;
    if(input>1){
        output = 32767;
    }else if(input<-1){
        output = -32768;
    }else
        output = 32768*input;
    return output;
}
//////////////
///////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*JNIEXPORT jboolean JNICALL
Java_com_google_sample_echo_MainActivity_createSLBufferQueueAudioPlayer(JNIEnv *env, jclass type) {
    SampleFormat sampleFormat;
    memset(&sampleFormat, 0, sizeof(sampleFormat));
    sampleFormat.pcmFormat_ = (uint16_t)engine.bitsPerSample_;
    sampleFormat.framesPerBuf_ = engine.fastPathFramesPerBuf_;

    // SampleFormat.representation_ = SL_ANDROID_PCM_REPRESENTATION_SIGNED_INT;
    sampleFormat.channels_ = (uint16_t)engine.sampleChannels_;
    sampleFormat.sampleRate_ = engine.fastPathSampleRate_;

    engine.player_ = new AudioPlayer(&sampleFormat, engine.slEngineItf_);
    assert(engine.player_);
    if(engine.player_ == nullptr)
        return JNI_FALSE;

    engine.player_->SetBufQueue(engine.recBufQueue_, engine.freeBufQueue_);
    engine.player_->RegisterCallback(EngineService, (void*)&engine);

    return JNI_TRUE;
}

JNIEXPORT void JNICALL
Java_com_google_sample_echo_MainActivity_deleteSLBufferQueueAudioPlayer(JNIEnv *env, jclass type) {
    if(engine.player_) {
        delete engine.player_;
        engine.player_= nullptr;
    }
}



JNIEXPORT void JNICALL
Java_com_google_sample_echo_MainActivity_deleteAudioRecorder(JNIEnv *env, jclass type) {
    if(engine.recorder_)
        delete engine.recorder_;

    engine.recorder_ = nullptr;
}

JNIEXPORT void JNICALL
Java_com_google_sample_echo_MainActivity_startPlay(JNIEnv *env, jclass type) {

    engine.frameCount_  = 0;
    /*
     * start player: make it into waitForData state
     */
/*   if(SL_BOOLEAN_FALSE == engine.player_->Start()){
       LOGE("====%s failed", __FUNCTION__);
       return;
   }
   engine.recorder_->Start();
}

JNIEXPORT void JNICALL
Java_com_google_sample_echo_MainActivity_stopPlay(JNIEnv *env, jclass type) {
   engine.recorder_->Stop();
   engine.player_ ->Stop();

   delete engine.recorder_;
   delete engine.player_;
   engine.recorder_ = NULL;
   engine.player_ = NULL;
}



uint32_t dbgEngineGetBufCount(void) {
   uint32_t count = engine.player_->dbgGetDevBufCount();
   count += engine.recorder_->dbgGetDevBufCount();
   count += engine.freeBufQueue_->size();
   count += engine.recBufQueue_->size();

   LOGE("Buf Disrtibutions: PlayerDev=%d, RecDev=%d, FreeQ=%d, "
                "RecQ=%d",
        engine.player_->dbgGetDevBufCount(),
        engine.recorder_->dbgGetDevBufCount(),
        engine.freeBufQueue_->size(),
        engine.recBufQueue_->size());
   if(count != engine.bufCount_) {
       LOGE("====Lost Bufs among the queue(supposed = %d, found = %d)",
            BUF_COUNT, count);
   }
   return count;
}

/*
* simple message passing for player/recorder to communicate with engine
*//*
bool EngineService(void* ctx, uint32_t msg, void* data ) {
    assert(ctx == &engine);
    switch (msg) {
        case ENGINE_SERVICE_MSG_KICKSTART_PLAYER:
            engine.player_->PlayAudioBuffers(PLAY_KICKSTART_BUFFER_COUNT);
            // we only allow it to call once, so tell caller do not call
            // anymore
            return false;
        case ENGINE_SERVICE_MSG_RETRIEVE_DUMP_BUFS:
            *(static_cast<uint32_t*>(data)) = dbgEngineGetBufCount();
            break;
        default:
            assert(false);
            return false;
    }

    return true;
}
*/



