package beast.evolution.datatype;

import beast.core.Description;
import beast.core.Input;
import beast.core.parameter.RealParameter;

@Description("Nucleotide error model from CellPhy paper with fixed parameters ADO and sequencing error")
public class NucleotideDiploid16WithError extends NucleotideDiploid16 implements DataTypeWithError {

    final public Input<RealParameter> deltaInput = new Input<>("delta", "the allelic dropout probability", Input.Validate.REQUIRED);
    final public Input<RealParameter> epsilonInput = new Input<>("epsilon", "the sequencing error probability",  Input.Validate.REQUIRED);

    private RealParameter delta;
    private RealParameter epsilon;

    protected double[][] errorMatrix;

    public NucleotideDiploid16WithError() {
        super();
    }

    @Override
    public void initAndValidate() {
        super.initAndValidate();
        delta = deltaInput.get();
        epsilon = epsilonInput.get();
        setupErrorMatrix();
    }

    @Override
    public void setupErrorMatrix() {
        setupErrorMatrix(delta.getValue(), epsilon.getValue());
    }

    private void setupErrorMatrix(double delta, double epsilon) {
        // true genotype is homozygous

        // true genotype is heterozygous

        double[][] matrix = {

        };
        errorMatrix = matrix;
    }

    public double getProbability(int observedState, int trueState) {
        if (isAmbiguousCode(observedState)) {
            return 1.0 / stateCount;
        } else {
            return errorMatrix[observedState][trueState];
        }
    }

    public double[] getProbabilities(int observedState) {
        double[] prob = new double[errorMatrix.length];
        for (int i = 0; i < errorMatrix.length; i++) {
            prob[i] = getProbability(observedState, i);
        }
        return prob;
    }

    @Override
    public String getTypeDescription() {
        return "nucleotideDiploid16WithError";
    }
}
