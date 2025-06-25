package phylonco.lphy.evolution.readcountmodel;

import lphy.core.logger.TextFileFormatted;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;

public class Integer2DMatrix implements TextFileFormatted{
    private Integer[][] matrix;

    public Integer2DMatrix() {}

    public Integer2DMatrix(Integer[][] matrix) {
        this.matrix = matrix;
    }

    public Integer[][] getMatrix() {
        return matrix;
    }

    public void setMatrix(Integer[][] matrix) {
        this.matrix = matrix;
    }

    public void setState (int taxon, int position, Integer state) {
        matrix[taxon][position] = state;
    }

    public Integer getState (int taxon, int position) {
        return matrix[taxon][position];
    }

    public Integer nTaxa () {return matrix.length;}

    public Integer nchar () {return matrix[0].length;}

    // not use old version at this moment, log new file instead
//    public String toString () {
//        StringBuilder sb = new StringBuilder();
//        sb.append("\"");
//        for (int i = 0; i < matrix.length; i++) {
//            for (int j = 0; j < matrix[i].length; j++) {
//                sb.append(matrix[i][j]);
//                if (j!=matrix[i].length - 1) {
//                    sb.append(",");
//                } else if (i != matrix.length - 1){
//                    sb.append("\"");
//                }
//            }
//        }
//        return sb.toString();
//    }

    @Override
    public List<String> getTextForFile() {
        // not using this at the moment
        return List.of();
    }

    @Override
    public void writeToFile(BufferedWriter writer) {
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[i].length; j++) {
                String line = String.valueOf(matrix[i][j]);
                if (j != matrix[i].length - 1) {
                    line = line + ",";
                } else if (i != matrix.length - 1) {
                    line = line + "\n";
                }
                try {
                    writer.write(line);
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }
        }
    }

    @Override
    public String getFileType() {
        return "_matrix.log";
    }
}
