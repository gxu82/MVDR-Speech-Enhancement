package com.example.abdullah.audiotwomics;

import android.media.AudioFormat;
import android.media.AudioManager;
import android.media.AudioTrack;
import android.util.Log;

import java.io.FileOutputStream;
import java.io.IOException;
import java.util.Arrays;
import java.util.concurrent.BlockingQueue;

/**
 * Created by Abdullah on 9/4/2017.
 */

public class AudioSaver implements Runnable {
    private FileOutputStream audioWriter;
    private BlockingQueue<AudioFrame> outputQueue;
    private boolean saveData = false;
    private AudioFrame outputFrame;
    private MainActivity main;
    private String outputName;
    private AudioTrack sound;
    private Thread writer;
    private static String TAG = AudioSaver.class.getSimpleName();

    public AudioSaver(String outputName, BlockingQueue<AudioFrame> outputQueue){
        main = Settings.getCallbackInterface();
        this.outputName = outputName;
        this.outputQueue = outputQueue;
        switch(Settings.debugLevel){
            case (0):{
                enablePCMoutput();
                break;
            }
            case (1):{
                enablePCMoutput();	// For Audio Save
                //Log.e(TAG, "onCreate initialize successfully!");
                break;
            }
        }

        if(Settings.playback){
            //int size = AudioTrack.getMinBufferSize(Settings.Fs, AudioFormat.CHANNEL_IN_STEREO, AudioFormat.ENCODING_PCM_16BIT);
            int size = Settings.frameSize*2;
            sound = new AudioTrack(AudioManager.STREAM_MUSIC, Settings.Fs, AudioFormat.CHANNEL_OUT_MONO, AudioFormat.ENCODING_PCM_16BIT, size, AudioTrack.MODE_STREAM);
            //sound = new AudioTrack(AudioManager.STREAM_MUSIC, Settings.Fs, AudioFormat.CHANNEL_OUT_STEREO, AudioFormat.ENCODING_PCM_16BIT, size, AudioTrack.MODE_STREAM);

        }
        Log.e(TAG, "onCreate initialize successfully!");
        writer = new Thread(this);
        writer.start();
    }

    private void enablePCMoutput() {
        try{
            if(main.getMode() == 2){
                audioWriter = new FileOutputStream(Utilities.getFile(outputName+ ".txt"));
                saveData = true;
            }
        }catch(IOException e){
            saveData = false;
            e.printStackTrace();
        }
    }

    public void run() {
        if(Settings.playback) {
            sound.play();
        }
        loop:while (true){
            outputFrame = null;
            try{
                outputFrame = outputQueue.take();
            }catch (InterruptedException e){
                Thread.currentThread().interrupt();
                e.printStackTrace();
            }
            //Log.e(TAG, "Audio takes successfully!");
            if(outputFrame!=null){
                if(outputFrame==Settings.STOP){
                    //if(outputFrame.getStopFlag()){
                    break loop;
                }
                if(Settings.playback) {
                    sound.write(outputFrame.getAudio(), 0, Settings.stepSize);
                }
                if(saveData){
                    try {
                        //audioWriter.write(Utilities.getByteArray(outputFrame.getAudio()));
                        audioWriter.write((Arrays.toString(outputFrame.getDebug()).replace("[", "").replace("]","")+", ").getBytes());
                        //audioWriter.write((Arrays.toString(outputFrame.getAudio())).getBytes());
                    } catch (IOException e) {
                        e.printStackTrace();
                    }
                }
                //Log.e(TAG, "Audio saves successfully!");
            }
        }
        main.done();
        if(sound != null){
            sound.pause();
            sound.flush();
            sound.release();
        }
        Log.e("AudioSaver","has been eliminated successfully!");
        if(saveData){
            try{
                audioWriter.flush();
                audioWriter.close();
            }catch (IOException e){
                e.printStackTrace();
            }
        }
    }
}
