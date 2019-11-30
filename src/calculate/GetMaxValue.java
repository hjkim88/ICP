/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package calculate;

/**
 *
 * @author KHJ
 */
public class GetMaxValue {

    private int[][] resultPair;
    private double[][] data;              // [gene pair1의 rowNum][gene pair2의 rowNum]
    private int th;                       // 상위 x개 Threshold
    private double[] dataPool;
    private double min;
    private int currentOrder;
    private int geneNum;
    private double[][][][][] meanResult;    // Inner-class clustering의 결과물    [gene1][gene2][class][cluster][axis]
    
    public GetMaxValue(double[][] d, int n, double[][][][][] meanResult) { //majority number
        data = d;
        th = n;
        this.meanResult = meanResult;
        initVariables();
        makeResultPair();
    }

    private void initVariables() {
        resultPair = new int [th][2];
        dataPool = new double[th];
        for(int i = 0; i < th; i++) {
            dataPool[i] = -1;
        }
        geneNum = data.length;
    }

    private void makeResultPair() {

        for(int i = 0; i < (geneNum-1); i++) {
            for(int j = (i+1); j < geneNum; j++) {
                getMin(dataPool);
                if(min < data[i][j]) {  //currentOrder = the index of minimum value in dataPool
                    dataPool[currentOrder] = data[i][j];
                    resultPair[currentOrder][0] = i;
                    resultPair[currentOrder][1] = j;
                }
                else if(min == data[i][j]) {
                    if(getSecondaryScore(meanResult[resultPair[currentOrder][0]][resultPair[currentOrder][1]]) < getSecondaryScore(meanResult[i][j])) {
                        dataPool[currentOrder] = data[i][j];
                        resultPair[currentOrder][0] = i;
                        resultPair[currentOrder][1] = j;
                    }
                }
            }
        }
    }

    private void getMin(double[] d) {
        int len = d.length;
        min = d[0];
        currentOrder = 0;
        
        for(int i = 1; i < len; i++) {
            if(d[i] == -1) {
                min = d[i];
                currentOrder = i;
                break;
            }
            else if(d[i] < min) {
                min = d[i];
                currentOrder = i;
            }
            else if(d[i] == min) {
                if(getSecondaryScore(meanResult[resultPair[i][0]][resultPair[i][1]]) < getSecondaryScore(meanResult[resultPair[currentOrder][0]][resultPair[currentOrder][1]])) {
                    min = d[i];
                    currentOrder = i;
                }
            }

        }
    }

    private double getSecondaryScore(double[][][] classifier) { // [class][cluster][axis]
        double x = 0, y = 0;
        double mx = 0, my = 0;
        int classNum = classifier.length;
        int clstrNum = classifier[0].length;

        for(int i = 0; i < classNum; i++) {
            for(int j = 0; j < clstrNum; j++) {
                x = x + classifier[i][j][0];
                y = y + classifier[i][j][1];
            }
        }

        mx = x / (double) (classNum * clstrNum);
        my = y / (double) (classNum * clstrNum);

        x = 0;
        y = 0;

        for(int i = 0; i < classNum; i++) {
            for(int j = 0; j < clstrNum; j++) {
                x = x + Math.pow(classifier[i][j][0] - mx, 2);
                y = y + Math.pow(classifier[i][j][1] - my, 2);
            }
        }
        
        return (x + y);
    }

    public int[][] getResultPair() {
        return resultPair;
    }

}
