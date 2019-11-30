/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package icp;

/**
 *
 * @author KHJ
 */
public class Main {

    int mostClass = 0;

    public Main() {
        initVariables();
    }

    private void initVariables() {

    }

    public void startICP() {
        double[][] ma1 = new file.LoadMA_Exl("D:/Data/prognosis/data2/GSE21034/S300/GSE21034_S300_ap19.xls").getMA();
        double[][] ma2 = new file.LoadMA_Exl("D:/Data/prognosis/data2/GSE21034/S300/GSE21034_S300_np131.xls").getMA();

        int c1Num = ma1[0].length;
        int c2Num = ma2[0].length;
        int sampleNum = c1Num + c2Num;
        int geneNum = ma1.length;

        if(c1Num > c2Num) {
            mostClass = 1;
        }
        else if(c1Num < c2Num) {
            mostClass = 2;
        }

        double[][][][][] meanResult;    // Inner-class clustering의 결과물    [gene1][gene2][class][cluster][axis]
        double[][] distanceScore = new double[geneNum][geneNum]; // 각 gene pair의 distance socre, 높을수록 좋음
        int[][] resultPair; // 상위 k개 gene pair
        int[][] result = new int[sampleNum][(2*global.Variables.candidateNum)+3];   // (0 : TP) (1 : FP) (2 : FN) (3 : TN)
        int[] candidatePrediction;

        for(int i = 0; i < global.Variables.clstrNum.length; i++) { // cluster number 만큼 돌림
        System.out.println("Cluster number = " + global.Variables.clstrNum[i]);
            for(int j = 0; j < sampleNum; j++) {    // sample number 만큼 돌림 - LOOCV 때문
                //System.out.println("Sample number = " + j + "/" + sampleNum);
                meanResult = new double[geneNum][geneNum][2][global.Variables.clstrNum[i]][2];
                for(int k = 0; k < (geneNum-1); k++) {
                    for(int l = (k+1); l < geneNum; l++) {  // k. l - gene pair 관련 연산
                        if(j < c1Num) { // 각 gene pair 마다 inner class clustering 의 결과물
                            meanResult[k][l][0] = new clustering.KMeans(global.Variables.clstrNum[i], makeTrainingData(ma1[k], j), makeTrainingData(ma1[l], j)).getMeanResult();
                            meanResult[k][l][1] = new clustering.KMeans(global.Variables.clstrNum[i], ma2[k], ma2[l]).getMeanResult();
                        }
                        else {
                            meanResult[k][l][0] = new clustering.KMeans(global.Variables.clstrNum[i], ma1[k], ma1[l]).getMeanResult();
                            meanResult[k][l][1] = new clustering.KMeans(global.Variables.clstrNum[i], makeTrainingData(ma2[k], j-c1Num), makeTrainingData(ma2[l], j-c1Num)).getMeanResult();
                        }
                        distanceScore[k][l] = new calculate.MeanDistance(meanResult[k][l]).getScore();  // 각 gene pair 마다 distance score 계산
                    }
                }

                resultPair = new calculate.GetMaxValue(distanceScore, global.Variables.candidateNum, meanResult).getResultPair();   // distance score를 바탕으로 각 loocv 마다 gene pair를 구함
                candidatePrediction = new int[global.Variables.candidateNum];   // loocv 마다 single prediction의 결과를 candidate number만큼 가지고 있음 - 1: c1, 2: c2

                for(int k = 0; k < global.Variables.candidateNum; k++) {    // 각 candidate(gene pair) 마다 single prediction을 하여 결과값 저장
                    if(j < c1Num) {
                        candidatePrediction[k] = singlePrediction(meanResult[resultPair[k][0]][resultPair[k][1]], ma1[resultPair[k][0]][j], ma1[resultPair[k][1]][j]);
                    }
                    else {
                        candidatePrediction[k] = singlePrediction(meanResult[resultPair[k][0]][resultPair[k][1]], ma2[resultPair[k][0]][j-c1Num], ma2[resultPair[k][1]][j-c1Num]);
                    }
                }

                for(int k = (-global.Variables.candidateNum-1); k <= (global.Variables.candidateNum+1); k++) {
                    if(j < c1Num) { // 각 LOOCV마다 최종 결과가 True인지 False 인지 계산 - (0 : TP) (1 : FP) (2 : FN) (3 : TN)
                        if(majorityVoting(candidatePrediction, 1, k) == true) {
                            result[j][k+global.Variables.candidateNum+1] = 0;
                        }
                        else {
                            result[j][k+global.Variables.candidateNum+1] = 2;
                        }
                    }
                    else {
                        if(majorityVoting(candidatePrediction, 2, k) == true) {
                            result[j][k+global.Variables.candidateNum+1] = 3;
                        }
                        else {
                            result[j][k+global.Variables.candidateNum+1] = 1;
                        }
                    }
                }

            }
            makeROC(result);
        }

    }

    public void getGenePair() {
        file.LoadMA_Exl loadMA = new file.LoadMA_Exl("D:/Data/prognosis/data2/GSE21034/S300/GSE21034_S300_ap19.xls");
        double[][] ma1 = loadMA.getMA();
        String[] geneSymbol = loadMA.getGeneSymbol();
        double[][] ma2 = new file.LoadMA_Exl("D:/Data/prognosis/data2/GSE21034/S300/GSE21034_S300_np131.xls").getMA();
        int geneNum = ma1.length;

        double[][][][][] meanResult;    // Inner-class clustering의 결과물    [gene1][gene2][class][cluster][axis]
        double[][] distanceScore = new double[geneNum][geneNum]; // 각 gene pair의 distance socre, 높을수록 좋음
        int[][] resultPair; // 상위 k개 gene pair

        for(int i = 0; i < global.Variables.clstrNum.length; i++) { // cluster number 만큼 돌림
            System.out.println("Cluster number = " + global.Variables.clstrNum[i]);
            meanResult = new double[geneNum][geneNum][2][global.Variables.clstrNum[i]][2];
            for(int j = 0; j < (geneNum-1); j++) {
                for(int k = (j+1); k < geneNum; k++) {  // k. l - gene pair 관련 연산
                    meanResult[j][k][0] = new clustering.KMeans(global.Variables.clstrNum[i], ma1[j], ma1[k]).getMeanResult();
                    meanResult[j][k][1] = new clustering.KMeans(global.Variables.clstrNum[i], ma2[j], ma2[k]).getMeanResult();
                    
                    distanceScore[j][k] = new calculate.MeanDistance(meanResult[j][k]).getScore();  // 각 gene pair 마다 distance score 계산
                }
            }

            resultPair = new calculate.GetMaxValue(distanceScore, global.Variables.candidateNum, meanResult).getResultPair();   // distance score를 바탕으로 각 loocv 마다 gene pair를 구함
            for(int j = 0; j < global.Variables.candidateNum; j++) {
                System.out.println(geneSymbol[resultPair[j][0]] + " " + geneSymbol[resultPair[j][1]]);
            }
        }

    }

    private void makeROC(int[][] result) {
        double AUC = 0;
        int sampleNum = result.length;
        int cutOffNum = result[0].length;

        double[] TP = new double[cutOffNum];
        double[] FP = new double[cutOffNum];
        double[] FN = new double[cutOffNum];
        double[] TN = new double[cutOffNum];

        double[] FPR = new double[cutOffNum];
        double[] TPR = new double[cutOffNum];

        for(int i = 0; i < cutOffNum; i++) {
            TP[i] = 0;
            FP[i] = 0;
            FN[i] = 0;
            TN[i] = 0;

            for(int j = 0; j < sampleNum; j++) {
                if(result[j][i] == 0) {
                    TP[i] ++;
                }
                else if(result[j][i] == 1) {
                    FP[i] ++;
                }
                else if(result[j][i] == 2) {
                    FN[i] ++;
                }
                else {
                    TN[i] ++;
                }
            }

            FPR[i] = FP[i] / (FP[i] + TN[i]);
            TPR[i] = TP[i] / (TP[i] + FN[i]);
        }

        System.out.println("Accuracy = " + ((TP[global.Variables.candidateNum+1] + TN[global.Variables.candidateNum+1]) / (double) sampleNum));

        AUC = calculateAUC(FPR, TPR);
        System.out.println("AUC = " + AUC);
        System.out.println("----------------------------------------------------");
        printCutOffPoints(FPR, TPR);
    }

    private void printCutOffPoints(double[] x, double[] y) {
        int cutOffNum = x.length;
        for(int i = 0; i < cutOffNum; i++) {
            System.out.println(x[i]);
        }
        System.out.println();
        for(int i = 0; i < cutOffNum; i++) {
            System.out.println(y[i]);
        }
    }

    private double calculateAUC(double[] x, double[] y) {
        double AUC = 0;
        int cutOffNum = x.length;

        for(int i = 0; i < (cutOffNum-1); i++) {
            AUC = AUC + ((y[i] + y[i+1]) * (x[i+1] - x[i]) / 2);
        }

        return AUC;
    }

    private boolean majorityVoting(int[] singleResult, int realClass, int plusVote) {
        boolean r = true;
        int candidateNum = singleResult.length;
        int cnt1 = plusVote, cnt2 = 0;
        int votedResult = 0;

        for(int i = 0; i < candidateNum; i++) {
            if(singleResult[i] == 1) {
                cnt1++;
            }
            else if(singleResult[i] == 2) {
                cnt2++;
            }
        }

        if(cnt1 > cnt2) {
            votedResult = 1;
        }
        else if(cnt1 < cnt2) {
            votedResult = 2;
        }
        else {
            votedResult = mostClass;
        }

        if(votedResult != realClass) {
            r = false;
        }

        return r;
    }

    private int singlePrediction(double[][][] classifier, double test_x, double test_y) {
        int r = 0;
        int clstrNum = classifier[0].length;
        double c1Min = getDistance(test_x, test_y, classifier[0][0][0], classifier[0][0][1]);
        double c2Min = getDistance(test_x, test_y, classifier[1][0][0], classifier[1][0][1]);
        double temp;

        for(int i = 1; i < clstrNum; i++) {
            temp = getDistance(test_x, test_y, classifier[0][i][0], classifier[0][i][1]);
            if(temp < c1Min) {
                c1Min = temp;
            }
            temp = getDistance(test_x, test_y, classifier[1][i][0], classifier[1][i][1]);
            if(temp < c2Min) {
                c2Min = temp;
            }
        }

        if(c1Min == c2Min) {
            c1Min = 0;
            c2Min = 0;
            for(int i = 0; i < clstrNum; i++) {
                c1Min = c1Min + getDistance(test_x, test_y, classifier[0][i][0], classifier[0][i][1]);
                c2Min = c2Min + getDistance(test_x, test_y, classifier[1][i][0], classifier[1][i][1]);
            }
        }

        if(c1Min < c2Min) {
            r = 1;
        }
        else if(c1Min > c2Min) {
            r = 2;
        }
        else {
            r = mostClass;
        }

        return r;
    }

    private double getDistance(double x1, double y1, double x2, double y2) {
        double d = 0;

        d = Math.pow((x1 - x2), 2) + Math.pow((y1 - y2), 2); //계산량을 줄이기 위해 sqrt는 하지 않음

        return d;
    }

    private double[] makeTrainingData(double[] d, int currentSample) {
        int len = d.length;
        double[] ma = new double[len-1];
        boolean isCurrent = false;

        for(int i = 0; i < len; i++) {
            if(i == currentSample) {
                isCurrent = true;
            }
            else if(isCurrent == false){
                ma[i] = d[i];
            }
            else {
                ma[i-1] = d[i];
            }
        }

        return ma;
    }

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        new icp.Main().startICP();
        //new icp.Main().getGenePair();
    }
}
