package com.example.abdullah.audiotwomics;

import android.Manifest;
import android.app.Activity;
import android.content.Context;
import android.content.pm.PackageManager;
import android.graphics.Color;
import android.os.Bundle;
import android.support.v4.app.ActivityCompat;
import android.text.method.ScrollingMovementMethod;
import android.widget.TextView;

import android.support.v7.app.AppCompatActivity;
import android.os.Bundle;
import android.util.Log;
import android.view.View;
import android.widget.Button;
import android.widget.EditText;
import android.widget.ImageButton;
import android.widget.RadioButton;
import android.widget.RadioGroup;
import android.widget.TextView;

import java.io.File;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.atomic.AtomicBoolean;

import static android.content.ContentValues.TAG;

public class MainActivity extends AppCompatActivity
        implements ActivityCompat.OnRequestPermissionsResultCallback {
    private int mode;

    private AtomicBoolean processingFlag;
    private RadioButton rbBabble;
    public TextView ake_text;

    private WaveRecorder waveRecorder;

    private boolean isEndFire=true;
    private boolean isBroadside=true;
    private boolean isOn=false;

    // private AudioReader audioReader;
    private BlockingQueue<AudioFrame> inputFrames;

    private BlockingQueue<AudioFrame> processedFrames;
    private static final int AUDIO_ECHO_REQUEST = 0;
    // Used to load the 'native-lib' library on application startup.
    static {
        System.loadLibrary("native-lib");
    }

    @Override
    protected void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);
        setContentView(R.layout.activity_main);

        Settings.setCallbackInterface(this);

        ake_text = (TextView)findViewById(R.id.textView11);
        ake_text.setMovementMethod(new ScrollingMovementMethod());
        inputFrames = new ArrayBlockingQueue<AudioFrame>(Settings.queueSize);
        processedFrames = new ArrayBlockingQueue<AudioFrame>(Settings.queueSize);
        processingFlag = new AtomicBoolean();
        Settings.ake_txt = ake_text;
        enableButtons(false);
        //enbaleSettings(false);
        setupListeners();
        // Example of a call to a native method
        //TextView tv = (TextView) findViewById(R.id.sample_text);
        //tv.setText(stringFromJNI());
    }

    private void enableButtons(boolean flag) {
        ((ImageButton)findViewById(R.id.ib_mic)).setEnabled(!flag);
        ((ImageButton)findViewById(R.id.ib_stop)).setEnabled(flag);
        ((Button)findViewById(R.id.ake_btn_endfire)).setEnabled(!flag);
        ((Button)findViewById(R.id.ake_btn_broadside)).setEnabled(!flag);
        ((Button)findViewById(R.id.ake_btn_onoff)).setEnabled(!flag);

    }

    private void setupListeners(){
        ((ImageButton)findViewById(R.id.ib_mic)).setOnClickListener(buttonClick);
        ((ImageButton)findViewById(R.id.ib_stop)).setOnClickListener(buttonClick);
        ((Button)findViewById(R.id.ake_btn_endfire)).setOnClickListener(buttonClick);
        ((Button)findViewById(R.id.ake_btn_broadside)).setOnClickListener(buttonClick);
        ((Button)findViewById(R.id.ake_btn_onoff)).setOnClickListener(buttonClick);

    }

    private View.OnClickListener buttonClick = new View.OnClickListener(){
        public void onClick(View v){
            switch (v.getId()){
                case R.id.ib_mic:{
                    //tv.setText(stringFromJNI());
                    StartRecord();

                    String output_cnt = Integer.toString((int)Settings.counter);
                    String output_timer = Float.toString(Settings.timer*1000);

                   // if ((int) Settings.counter % 10==0){
                        ake_text.setTextColor(Color.RED);
                        ake_text.append("\n"+"FFT calculation of frame "+ output_cnt + " takes "+output_timer+ "ms"  );
                    //}

                    break;
                }
                case R.id.ib_stop:{
                    StopRecord();
                    break;
                }
                case R.id.ake_btn_endfire:{
                    if(isEndFire){
                        ((Button) findViewById(R.id.ake_btn_endfire)).setBackgroundColor(Color.RED);
                        isEndFire=false;
                        ((Button) findViewById(R.id.ake_btn_broadside)).setBackgroundColor(Color.GREEN);
                        isBroadside=true;
                    }else {
                        ((Button) findViewById(R.id.ake_btn_endfire)).setBackgroundColor(Color.GREEN);
                        isEndFire=true;
                        ((Button) findViewById(R.id.ake_btn_broadside)).setBackgroundColor(Color.RED);
                        isBroadside=false;
                    }
                    changefilterdirection();
                    break;
                }
                case R.id.ake_btn_broadside:{
                    if(isBroadside){
                        ((Button) findViewById(R.id.ake_btn_endfire)).setBackgroundColor(Color.GREEN);
                        isEndFire=true;
                        ((Button) findViewById(R.id.ake_btn_broadside)).setBackgroundColor(Color.RED);
                        isBroadside=false;
                    }else {
                        ((Button) findViewById(R.id.ake_btn_endfire)).setBackgroundColor(Color.RED);
                        isEndFire=false;
                        ((Button) findViewById(R.id.ake_btn_broadside)).setBackgroundColor(Color.GREEN);
                        isBroadside=true;
                    }
                    changefilterdirection();
                    break;
                }
                case R.id.ake_btn_onoff:{
                    if(isOn){
                        //((Button) findViewById(R.id.ake_btn_endfire)).setBackgroundColor(Color.GREEN);
                        isOn=false;
                        ((Button) findViewById(R.id.ake_btn_onoff)).setText("OFF");
                    }else {
                      //  ((Button) findViewById(R.id.ake_btn_endfire)).setBackgroundColor(Color.RED);
                        isOn=true;
                        ((Button) findViewById(R.id.ake_btn_onoff)).setText("ON");
                    }
                    activatefiltering();
                    break;
                }
            }
        }
    };

    private void StopRecord() {
        switch(mode){
            case 1:{
                waveRecorder.stopRecording();
                waveRecorder = null;
               // Settings.ake_txt.setTextColor(Color.RED);

                Settings.ake_txt.append("\n"+"It is stopped succesfully"  );
                break;
            }
        }
    }

    private void StartRecord() {
        if (ActivityCompat.checkSelfPermission(this, Manifest.permission.RECORD_AUDIO) !=
                PackageManager.PERMISSION_GRANTED) {

            ActivityCompat.requestPermissions(
                    this,
                    new String[] { Manifest.permission.RECORD_AUDIO },
                    AUDIO_ECHO_REQUEST);
            return;
        }

        if (ActivityCompat.checkSelfPermission(this, Manifest.permission.WRITE_EXTERNAL_STORAGE) !=
                PackageManager.PERMISSION_GRANTED) {

            ActivityCompat.requestPermissions(
                    this,
                    new String[] { Manifest.permission.WRITE_EXTERNAL_STORAGE },
                    AUDIO_ECHO_REQUEST);
            return;
        }

        if (ActivityCompat.checkSelfPermission(this, Manifest.permission.MODIFY_AUDIO_SETTINGS) !=
                PackageManager.PERMISSION_GRANTED) {

            ActivityCompat.requestPermissions(
                    this,
                    new String[] { Manifest.permission.MODIFY_AUDIO_SETTINGS },
                    AUDIO_ECHO_REQUEST);
            return;
        }

        inputFrames.clear();

        processedFrames.clear();

        new Processing(inputFrames,processedFrames);


        waveRecorder = new WaveRecorder(inputFrames);

        if(Settings.debugLevel!=2){
            mode = 1;
            new AudioSaver(getString(R.string.appDirectory)+ File.separator+System.currentTimeMillis(),processedFrames);
        }

        //startPlay();

        processingFlag.set(true);
        enableButtons(true);
    }

    public int getMode() {
        return mode;
    }

    public void done() {
        if(processingFlag.getAndSet(false)) {
            runOnUiThread(
                    new Runnable() {
                        public void run() {
                            enableButtons(false);
                        }
                    }
            );
        }
    }

    private void changefilterdirection() {


        Settings.changefilterdirection(isEndFire,isBroadside);

    }

    private void activatefiltering() {


        Settings.activatefiltering(isOn);

    }
    public void updateGUI()
    {
        String output_cnt = Integer.toString((int)Settings.counter);
        String output_timer = Float.toString(Settings.timer);
        Settings.ake_txt.setTextColor(Color.RED);
         if ((int) Settings.counter % 10==0){

        //ake_text.setText("\n"+"FFT calculation of frame "+ output_cnt + "takes"+output_timer+ "ms"  );
        Settings.ake_txt.append("\n"+"FFT calculation of frame "+ output_cnt + " takes "+output_timer+ " ms"  );
        }
    }
    /**
     * A native method that is implemented by the 'native-lib' native library,
     * which is packaged with this application.
     */
  //  public native String stringFromJNI();
}
