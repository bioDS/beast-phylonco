package phylonco.lphy.evolution.readcountmodel;

public class Integer2DMatrix {
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

}
