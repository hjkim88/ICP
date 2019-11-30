/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package calculate;

/**
 *
 * @author KHJ
 */
public class MeanDistance {
    
    private double[][][] meanResult;    //  [class][cluster][axis]
    private int clstrNum;
    private double distance;

    public MeanDistance(double[][][] d) {  // [class][cluster][axis]
        meanResult = d;
        initVariables();
        calMeanDistance();
    }

    private void initVariables() {
        clstrNum = meanResult[0].length;
    }

    private double getDistance(double x1, double y1, double x2, double y2) {
        double d = 0;

        d = Math.pow((x1 - x2), 2) + Math.pow((y1 - y2), 2); //계산량을 줄이기 위해 sqrt는 하지 않음

        return d;
    }

    private void calMeanDistance() {
        distance = 0;

        for(int i = 0; i < clstrNum; i++) {
            for(int j = 0; j < clstrNum; j++) {
                distance = distance + getDistance(meanResult[0][i][0], meanResult[0][i][1], meanResult[1][j][0], meanResult[1][j][1]);
            }
        }
    }

    public double getScore() {
        return distance;
    }
}
