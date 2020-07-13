package com.example.abdullah.audiotwomics;

import android.os.Environment;

import java.io.File;

/**
 * Created by Abdullah on 9/4/2017.
 */

public class Utilities {
    public static File getFile(String outName){
        File outfile = new File(Environment.getExternalStorageDirectory().getPath(), outName);
        if (!outfile.exists()) {
            try {
                outfile.createNewFile();
            } catch (Exception e) {}
        }
        return outfile;
    }

    public static File prepareDirectory(String outName){
        File outDir = new File(Environment.getExternalStorageDirectory().getPath(), outName);
        if (outDir.isFile()) {
            if(outDir.delete()){
                outDir.mkdirs();
            }
        }
        if(!outDir.isDirectory()){
            outDir.mkdirs();
        }
        return outDir;
    }
    public static byte[] getByteArray(short[] samples) {

        byte[] bytes = new byte[samples.length*2];

        for (int i=0;i<samples.length;i++) {
            bytes[2*i] = (byte)(samples[i] & 0x00FF);
            bytes[(2*i)+1] = (byte)(samples[i]>>8);
        }
        return bytes;
    }
}