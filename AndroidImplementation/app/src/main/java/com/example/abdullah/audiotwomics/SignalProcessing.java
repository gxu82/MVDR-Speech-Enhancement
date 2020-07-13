package com.example.abdullah.audiotwomics;

/**
 * Created by Abdullah on 9/4/2017.
 */

public class SignalProcessing {

    private long pointer;
    public SignalProcessing(){
        pointer = paramInitilization(Settings.Fs,Settings.stepSize,Settings.frameSize, Settings.stepFactor, Settings.noisetype,Settings.isEndFire,Settings.isOn);
    }

    public void frameProcess(short[] input){
        realtimeProcessing(pointer,input);
    }

    public float getTime(){
        return getTime(pointer);
    }

    public float getComputeTime(){
        return getComputeTime(pointer);
    }

    public float getFilteringTime(){
        return getFilteringTime(pointer);
    }

    public void release(){
        paramElimination(pointer);
    }

    public short[] getSoundOutput(){
        return soundOutput(pointer,Settings.debugLevel);
    }

    public float[] getTiming(){
        return FFTtiming(pointer,Settings.debugLevel);
    }

    public short[] getDataOutput(){
        return dataOutput(pointer,Settings.debugLevel);
    }


    public static native long paramInitilization(int frequency, int stepSize,int frameSize,float stepFactor, int noisetype, boolean isEnchanced, boolean isOn);
    public static native void paramElimination(long memoryPointer);
    public static native float getTime(long memoryPointer);
    public static native float getComputeTime(long memoryPointer);
    public static native float getFilteringTime(long memoryPointer);
    public static native void realtimeProcessing(long memoryPointer,short[] input);
    public static native short[] soundOutput(long memoryPointer,int outputSelection);

    public static native float[] FFTtiming(long memoryPointer,int outputSelection);

    public static native short[] dataOutput(long memoryPointer,int outputSelection);
    //public static native void startPlay();


}
