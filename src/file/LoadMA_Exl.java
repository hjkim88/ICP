/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package file;

import java.io.File;
import java.io.IOException;
import jxl.Workbook;
import jxl.Sheet;
import jxl.NumberCell;

/**
 *
 * @author KHJ
 */
public class LoadMA_Exl {

    private double[][] ma;
    private String[] geneSymbol;

    public LoadMA_Exl(String maPath) {
        initVariables();
        read(maPath);
    }

    private void initVariables() {
        if(global.Variables.debuggingMode == true) {
            System.out.println("file.LoadMA_Exl() ... Begins.");
        }
    }

    private void read(String maPath) {
        try {
            Workbook workbook = Workbook.getWorkbook(new File(maPath));
            Sheet sheet = workbook.getSheet(0);

            int rowNum = sheet.getRows()-2;
            int colNum = sheet.getColumns()-1;

            ma = new double[rowNum][colNum];
            geneSymbol = new String[rowNum];
            NumberCell nc;

            for(int i = 0; i < rowNum; i++) {
                geneSymbol[i] = sheet.getCell(0, i+2).getContents();
                for(int j = 0; j < colNum; j++) {
                    nc = (NumberCell) sheet.getCell(j+1, i+2);
                    nc.getNumberFormat().setMaximumFractionDigits(10);
                    ma[i][j] = nc.getValue();
                }
            }

            workbook.close();
        }
        catch(IOException ioe) {
            if(global.Variables.debuggingMode == true) {
                System.out.println("file.LoadMA_Exl().read().IOException");
            }
        }
        catch(jxl.JXLException jxle) {
            if(global.Variables.debuggingMode == true) {
                System.out.println("file.LoadMA_Exl().read().JXLException");
            }
        }
    }

    public String[] getGeneSymbol() {
        if(global.Variables.debuggingMode == true) {
            System.out.println("file.LoadMA_Exl().getGeneSymbol()");
        }
        return geneSymbol;
    }
    
    public double[][] getMA() {
        if(global.Variables.debuggingMode == true) {
            System.out.println("file.LoadMA_Exl().getMA()");
        }
        return ma;
    }

}
