package com.example.abdullah.audiotwomics;

/**
 * Created by Abdullah Kucuk on 9/4/2017.
 */

import android.media.AudioFormat;
import android.media.AudioManager;
import android.media.AudioRecord;
import android.media.AudioTrack;
import android.media.MediaRecorder;
import android.util.Log;

import java.util.concurrent.BlockingQueue;
import java.util.concurrent.atomic.AtomicBoolean;

public class WaveRecorder implements Runnable {
    private short[] abuf_short;
    private AtomicBoolean recordFlag;
    private int recBufLen;
    private AudioRecord recorder;
    //private AudioTrack player;
    private BlockingQueue<AudioFrame> output;
    private Thread recThread;
    private static String TAG = WaveRecorder.class.getSimpleName();
    public WaveRecorder(BlockingQueue<AudioFrame> output){
        this.output = output;
        recordFlag = new AtomicBoolean();
        // recBufLen = AudioRecord.getMinBufferSize(Settings.Fs, AudioFormat.CHANNEL_IN_STEREO, AudioFormat.ENCODING_PCM_16BIT);
        recBufLen=Settings.stepSize*2;
        recorder = new AudioRecord(MediaRecorder.AudioSource.MIC,Settings.Fs, AudioFormat.CHANNEL_IN_STEREO, AudioFormat.ENCODING_PCM_16BIT,recBufLen);
       // player = new AudioTrack(AudioManager.STREAM_MUSIC, Settings.Fs, AudioFormat.CHANNEL_OUT_STEREO,AudioFormat.ENCODING_PCM_16BIT, recBufLen*2, AudioTrack.MODE_STREAM);
        Log.e(TAG, "onCreate initialize successfully!");
        recThread = new Thread(this);
        recThread.start();
    }

    public void stopRecording(){
        recordFlag.set(false);
    }

    public boolean recordFlag(){
        return recordFlag.get();
    }

    public void run() {
        //android.os.Process.setThreadPriority(android.os.Process.THREAD_PRIORITY_URGENT_AUDIO);
        try{
            if(recorder!=null){
                recorder.stop();
                recorder.startRecording();
             //   if(player.getPlayState() != AudioTrack.PLAYSTATE_PLAYING)
             //      player.play();
                recordFlag.set(true);

            }
        }catch(IllegalStateException e){
            recordFlag.set(false);
        }
        try{
            loop:while(true){
                if(recordFlag.get()){
                    abuf_short = new short[2*Settings.stepSize];
                    recorder.read(abuf_short,0,2*Settings.stepSize);
                 //   player.write(abuf_short,0,2*Settings.stepSize);
                    output.put(new AudioFrame(abuf_short,false));
                    //Log.e(TAG, "Recording successfully!");
                }else{
                    output.put(Settings.STOP);
                    break loop;
                }
            }
        }catch (InterruptedException e){
            Thread.currentThread().interrupt();
            e.printStackTrace();
        }

        try{
            if(recorder!=null){
                recorder.stop();
                recorder.release();
                recorder=null;
               // player.stop();
               // player.release();
               // player=null;
                Log.e("WaveRecorder"," has been eliminated successfully!");
            }
        }catch (IllegalStateException e){
            e.printStackTrace();
        }
    }
}
