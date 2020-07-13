package com.example.abdullah.audiotwomics;

import android.app.Activity;
import android.util.Log;

import java.util.concurrent.BlockingQueue;

/**
 * Created by Abdullah Kucuk on 9/4/2017.
 */

public class Processing extends Activity implements Runnable {




    MainActivity myActivity = new MainActivity();


    private BlockingQueue<AudioFrame> input;

    private BlockingQueue<AudioFrame> output;
    private SignalProcessing signalProcessor;
    private Thread speechThread;
    //private int counter;
    private float[] temp_angle;
    private short[] temp_output;
    private float counter;
    private float timer;
    private static String TAG = Processing.class.getSimpleName();
    public Processing(BlockingQueue<AudioFrame> input, BlockingQueue<AudioFrame> output){
        this.input = input;

        this.output = output;
        signalProcessor = new SignalProcessing();
        speechThread = new Thread(this);
        speechThread.start();
        Log.e(TAG, "onCreate initialize successfully!");
    }
    @Override
    public void run() {
        try{
            loop:while(true){
                AudioFrame currentFrame = null;



                currentFrame = input.take();



                if(currentFrame==Settings.STOP){
                   // Log.e(TAG,"Stage 1 compute time:"+signalProcessor.getFilteringTime()+" ms");
                   // Log.e(TAG,"Stage 2 compute time:"+signalProcessor.getComputeTime()+" ms");
                    signalProcessor.release();
                    output.put(currentFrame);
                    break loop;
                }
                //Log.e(TAG, "frame taken successfully!");

                signalProcessor.frameProcess(currentFrame.getAudio());
                //Log.e(TAG, "frame Processing successfully!");

                if(Settings.playback){
                    currentFrame.setAduio(signalProcessor.getSoundOutput());
                    temp_angle = signalProcessor.getTiming();
                    counter = temp_angle[0];
                    Settings.counter = counter;
                    timer = temp_angle[1];
                    Settings.timer = timer;
                    //angle = temp_angle[0];

                    myActivity.runOnUiThread(new Runnable() {
                        public void run() {
                            myActivity.updateGUI();
                        }
                    });



                }
                /*
                if(Settings.debugLevel==0){
                    currentFrame.setDebug(signalProcessor.getDataOutput());
                }
                */
                output.put(currentFrame);
            }
        }catch(InterruptedException e){
            Thread.currentThread().interrupt();
            e.printStackTrace();
        }
    }
}