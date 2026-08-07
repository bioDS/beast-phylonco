package phylonco.beast.evolution.errormodel;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.evolution.datatype.DataType;
import beast.base.spec.type.RealScalar;
import phylonco.beast.evolution.datatype.NucleotideDiploid10;

import static beast.base.evolution.datatype.DataType.GAP_CHAR;
import static beast.base.evolution.datatype.DataType.MISSING_CHAR;

/**
 * Implements the GT10 error model for diploid genotypes from Kozlov et al. (2021)
 *
 * CellPhy: accurate and fast probabilistic inference of single-cell phylogenies from scDNA-seq data
 * https://doi.org/10.1101/2020.07.31.230292
 *
 */
@Description("GT10 diploid error model from CellPhy paper")
public class GT10ErrorModel extends ErrorModel {

    final public Input<RealScalar> deltaInput = new Input<>("delta", "the allelic dropout probability", Input.Validate.REQUIRED);
    final public Input<RealScalar> epsilonInput = new Input<>("epsilon", "the sequencing error probability",  Input.Validate.REQUIRED);

    private RealScalar delta;
    private RealScalar epsilon;

    @Override
    public void initAndValidate() {
        super.initAndValidate();

        delta = deltaInput.get();
        epsilon = epsilonInput.get();

        if (!(datatype instanceof NucleotideDiploid10)) {
            throw new IllegalArgumentException("datatype must be of class NucleotideDiploid10");
        }

        if (updateMatrix) {
            setupErrorMatrix();
            updateMatrix = false;
        }
    }

    @Override
    public double getProbability(int observedState, int trueState) {
        double prob;
        String observedStr = datatype.getCharacter(observedState);
        String gap = Character.toString(GAP_CHAR);
        String missing = Character.toString(MISSING_CHAR);
        if (observedStr.equals(gap) || observedStr.equals(missing)) {
            // gap or missing code
            prob = 1.0;
        } else if (datatype.isAmbiguousCode(observedState)) {
            // ambiguous code for more than one state (not gap or missing)
            int[] codes = datatype.getStatesForCode(observedState);
            prob = 0.0;
            for (int i = 0; i < codes.length; i++) {
                prob += getProbabilityUnambiguous(codes[i], trueState);
            }
        } else {
            // unambiguous code
            prob = getProbabilityUnambiguous(observedState, trueState);
        }

        if (prob > 1.0) {
            throw new RuntimeException("Error in GT16ErrorModel: Tip partial likelihood cannot exceed 1.0!");
        }

        return prob;
    }

    // the error probabilities are the same as GT16 except for scenario II and VI which are updated to
    // scenario II: P (ab | aa) = (1 – delta) * (1/3) * epsilon
    // scenario VI: P(ac | ab) = P (ca | ab) = P(cb | ab) = P(bc | ab)
    private double getProbabilityUnambiguous(int observedState, int trueState) {
        double d = delta.get();
        double e = epsilon.get();

        NucleotideDiploid10 gt10Datatype = (NucleotideDiploid10)(datatype);
        String trueGenotype = gt10Datatype.getGenotype(trueState);
        String observedGenotype = gt10Datatype.getGenotype(observedState);

        char trueFirst = trueGenotype.charAt(0);
        char trueSecond = trueGenotype.charAt(1);
        char observedFirst = observedGenotype.charAt(0);
        char observedSecond = observedGenotype.charAt(1);

        double prob = 0.0;

        if (trueFirst == trueSecond) {
            // true state homozygous
            if (observedState == trueState) {
                // scenario I
                // P(aa | aa) = (1 - epsilon) + (1/2) * delta * epsilon
                prob = 1 - e + (1.0 / 2.0) * d * e;
            } else if (observedFirst == trueFirst || observedSecond == trueSecond) {
                // scenario II updated to
                // P(ab | aa) =
                // P(ba | aa) = (1 - delta) * (1/3) * epsilon
                prob = (1 - d) * (1.0 / 3.0) * e;
            } else if (observedFirst == observedSecond) {
                // scenario III
                // P(bb | aa) = (1/6) * delta * epsilon
                prob = (1.0 / 6.0) * d * e;
            } else {
                // P(bc | aa) = 0
                prob = 0.0;
            }
        } else {
            // true state heterozygous
            if (observedState == trueState) {
                // scenario VII
                // P(ab | ab) = (1 - delta) * (1 - epsilon)
                prob = (1 - d) * (1 - e);
            } else if (observedFirst == observedSecond) {
                // observed is homozygous
                if (observedFirst == trueFirst || observedSecond == trueSecond) {
                    // scenario IV
                    // P(aa | ab) =
                    // P(bb | ab) = (1/2) * delta + (1/6) * epsilon - (1/3) * delta * epsilon
                    prob = (1.0 / 2.0) * d + (1.0 / 6.0) * e - (1.0 / 3.0) * d * e;
                } else {
                    // scenario V
                    // P(cc | ab) = (1/6) * delta * epsilon
                    prob = (1.0 / 6.0) * d * e;
                }
            } else {
                // observed is heterozygous
                if (observedFirst == trueFirst || observedSecond == trueFirst ||
                        observedSecond == trueSecond || observedFirst == trueSecond  ) {
                    // scenario VI
                    // P(ac | ab) = P(ca | ab) = P(cb | ab) = P(bc | ab)
                    // = (1 - delta) * (1/6) * epsilon
                    prob = (1 - d) * (1.0 / 6.0) * e;
                } else {
                    // P(cd | ab) = 0
                    prob = 0.0;
                }
            }
        }

        return prob;
    }

    @Override
    public double[] getProbabilities(int observedState) {
        if (updateMatrix) {
            setupErrorMatrix();
            updateMatrix = false;
        }
        return errorMatrix[observedState];
    }

    @Override
    public void setupErrorMatrix() {
        if (errorMatrix == null) {
            errorMatrix = new double[datatype.mapCodeToStateSet.length][datatype.getStateCount()];
        }
        for (int trueState = 0; trueState < datatype.getStateCount(); trueState++) {
            for (int observedState = 0; observedState < datatype.mapCodeToStateSet.length; observedState++) {
                // rows are observed states X, columns are true states Y
                errorMatrix[observedState][trueState] = getProbability(observedState, trueState);
            }
        }
    }

    @Override
    public boolean canHandleDataType(DataType datatype) {
        return datatype instanceof NucleotideDiploid10;
    }

}
