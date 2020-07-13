package com.example.abdullah.audiotwomics;

import android.widget.EditText;
import android.widget.TextView;

/**
 * Created by Abdullah on 9/4/2017.
 */

public class Settings {
    public static int Fs = 48000;
    public static float stepTime = 20.0f;
    public static float frameTime = 40.0f;
    public static int stepSize = Math.round(Fs*stepTime*0.001f);
    public static int frameSize = Math.round(Fs*frameTime*0.001f);
    public static float stepFactor = 0.5f;
    public static final int queueSize = 10;


    public static boolean isEndFire = true;
    public static boolean isBroadSide = false;
    public static boolean isOn = false;

    public static boolean playback = true;
    public static int noisetype = 1;
    public static int debugLevel = 0;

    public static TextView ake_txt;
    public static EditText ake_timer_txt;


    public static float counter = 0;
    public static float timer = 0;


    private static MainActivity main;
    public static final AudioFrame STOP = new AudioFrame(new short[] {1,2,4,8},true);
    public static void setCallbackInterface(MainActivity uiInterface) {
        main = uiInterface;
    }
    public static MainActivity getCallbackInterface() {
        return main;
    }

    public static boolean setupStepSize(float new_stepsize){
        if (new_stepsize>0){
            stepFactor = new_stepsize;
        }
        return false;
    }



    public static boolean setupFrameTime(float new_frametime){
        if (new_frametime>0){
            frameTime = new_frametime;
            stepTime = frameTime/2;
            frameSize = Math.round(Fs*frameTime*0.001f);
            stepSize = Math.round(Fs*stepTime*0.001f);
            return true;
        }
        return false;
    }

    public static boolean setupNoiseType(int new_noisetype){
        if (new_noisetype!=noisetype){
            noisetype = new_noisetype;
            return true;
        }
        return false;
    }

    public static boolean setupPlayback(int new_playback){
        if(debugLevel!=new_playback){
            debugLevel = new_playback;
            return true;
        }
        return false;
    }
    public static void changefilterdirection(boolean new_isEndFire,boolean new_isBroadside){

        isEndFire = new_isEndFire;
        isBroadSide = new_isBroadside;

        //return isEndFire;
    }
    public static void activatefiltering(boolean new_isOn){

        isOn = new_isOn;


        //return isEndFire;
    }

}
