/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package clustering;

import java.util.Random;

/**
 *
 * @author KHJ
 */
public class KMeans {

    private int k = 0;
    private double[] dataX;
    private double[] dataY;
    private double[] meanX;
    private double[] meanY;
    private int[] initKPool;
    private boolean isExist;
    private int[] cluster;
    private int[] oldCluster;
    private boolean isChanged;
    private double[][] meanResult;      // [cluster][axis]
    private int threshold;

    public KMeans(int clstrNum, double[] data1, double[] data2) {
        k = clstrNum;
        dataX = data1;
        dataY = data2;
        initVariables();
        startKMeans();
    }

    private void initVariables() {
        meanX = new double[k];
        meanY = new double[k];
        initKPool = new int[k];
        cluster = new int[dataX.length];
        oldCluster = new int[dataX.length];
        isExist = false;
        isChanged = true;
        meanResult = new double[k][2];
        threshold = global.Variables.kMeansThreshold;               // k-Means 반복하는 쓰레숄드 지정.
    }
    
    public double[][] getMeanResult() {
        return meanResult;
    }

    private void startKMeans() {
        //selectInitK_Random();
        selectInitK_Determined();
        
        int cnt = 0;
        while((isChanged == true) && (cnt < threshold)) {
            calculateDistanceOfAll();
            compareTwoCluster();
            if(isChanged == true) {
                clusterRebuilder();
                cnt++;
            }
        }
        for(int i = 0; i < k; i++) {
            meanResult[i][0] = meanX[i];
            meanResult[i][1] = meanY[i];
        }
    }

    private void selectInitK_Random() {
        for(int i = 0; i < k; i++) {
            initKPool[i] = (-1);
        }

        int n = 0;
        Random rand = new Random(System.currentTimeMillis());
        while(n < k) {
            int z = Math.abs(rand.nextInt(dataX.length));
            testIsExist(z);
            if(isExist == false) {
                initKPool[n] = z;
                n++;
            }
        }

        for(int j = 0; j < k; j++) {
            meanX[j] = dataX[initKPool[j]];
            meanY[j] = dataY[initKPool[j]];
        }
        
        for(int t = 0; t < cluster.length; t++) {
            cluster[t] = (-1);
        }

    }

    private void selectInitK_Determined() {
        int separator = (dataX.length - 1) / (k - 1);

        for(int i = 0; i < k; i++) {
            meanX[i] = dataX[i * separator];
            meanY[i] = dataY[i * separator];
        }

        for(int t = 0; t < cluster.length; t++) {
            cluster[t] = (-1);
        }
    }

    private void calculateDistanceOfAll() {
        for(int a = 0; a < cluster.length; a++) {
            oldCluster[a] = cluster[a];
        }
        double[] d = new double[k];
        double min = 0;
        for(int i = 0; i < dataX.length; i++) {
            for(int j = 0; j < k; j++) {
                d[j] = getDistance(dataX[i], dataY[i], meanX[j], meanY[j]);
                if(j == 0) {
                    min = d[j];
                }
                min = Math.min(min, d[j]);
            }
            for(int j = 0; j < k; j++) {
                if(d[j] == min) {
                    cluster[i] = j;
                }
            }
        }
    }

    private void clusterRebuilder() {
        double[] sumX = new double[k];
        double[] sumY = new double[k];
        int cnt[] = new int[k];

        for(int n = 0; n < k; n++) {
            cnt[n] = 0;
        }
        
        for(int i = 0; i < dataX.length; i++) {
            for(int j = 0; j < k; j++) {
                if(cluster[i] == j) {
                    sumX[j] = sumX[j] + dataX[i];
                    sumY[j] = sumY[j] + dataY[i];
                    cnt[j]++;
                }
            }
        }

        for(int t = 0; t < k; t++) {
            meanX[t] = (double) sumX[t] / cnt[t];
            meanY[t] = (double) sumY[t] / cnt[t];
        }
    }

    private void compareTwoCluster() {
        isChanged = false;
        for(int i = 0; i < cluster.length; i++) {
            if(cluster[i] != oldCluster[i]) {
                isChanged = true;
            }
        }
    }

    private double getDistance(double x1, double y1, double x2, double y2) {
        double d = 0;

        d = Math.pow((x1 - x2), 2) + Math.pow((y1 - y2), 2); // 계산량을 줄이기 위해 sqrt는 하지 않음

        return d;
    }

    private void testIsExist(int z) {
        isExist = false;
        for(int i = 0; i < initKPool.length; i++) {
            if(initKPool[i] == z) {
                isExist = true;
            }
        }
    }

}
