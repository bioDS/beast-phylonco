package phylonco.beast.evolution.substitutionmodel;

import beast.core.Description;
import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.evolution.datatype.DataType;
import beast.evolution.datatype.NucleotideMethylation;
import beast.evolution.substitutionmodel.ComplexColtEigenSystem;
import beast.evolution.substitutionmodel.ComplexSubstitutionModel;
import beast.evolution.substitutionmodel.SubstitutionModel;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.SingularValueDecomposition;

@Description("Covarion model for methylation data based on HKY nucleotide substitution model.")
public class MethylationHKY extends ComplexSubstitutionModel implements SubstitutionModel {
    public Input<RealParameter> kappaInput = new Input<RealParameter>(
            "kappa", "kappa parameter of the HKY model", Input.Validate.REQUIRED);
    public Input<RealParameter> alphaInput = new Input<RealParameter>(
            "alpha", "rate of methylation (C->MetC)", Input.Validate.REQUIRED);
    public Input<RealParameter> betaInput = new Input<RealParameter>(
            "beta", "rate of demethylation (MetC->C)", Input.Validate.REQUIRED);
    public Input<RealParameter> gammaInput = new Input<RealParameter>(
            "gamma", "rate of demethylation of MetC into Thymine", Input.Validate.REQUIRED);


    private RealParameter kappaPar;
    private RealParameter alphaPar;
    private RealParameter betaPar;
    private RealParameter gammaPar;

    public MethylationHKY() {
        ratesInput.setRule(Input.Validate.OPTIONAL);
        frequenciesInput.setRule(Input.Validate.OPTIONAL);
    }

    @Override
    public void initAndValidate() {
        nrOfStates = 6; // optionally nrOfStates = frequencies.getFreqs().length;
        kappaPar = kappaInput.get();
        alphaPar = alphaInput.get();
        betaPar = betaInput.get();
        gammaPar = gammaInput.get();

        rateMatrix = new double[nrOfStates][nrOfStates];

        // eigenSystem for eigen decomposition to get transition matrix in GeneralSubstitutionModel
        // couldn't set default value to be RobustEigenSystem otherwise
        //eigenSystem = new RobustEigenSystem(getStateCount());
        eigenSystem = new ComplexColtEigenSystem(nrOfStates, false, 1000000, 1000000);
    }

    @Override
    public boolean canHandleDataType(DataType dataType) {
        if (dataType instanceof NucleotideMethylation) {
            return true;
        } else {
            return false;
        }
    }



    @Override
    protected void setupRelativeRates() {
    }


    /** Get frequencies of stationary distribution
     *
     * Frequencies of stationary distribution are derived from observed frequencies of 0, 1 and W in data.
     * It is also assumed that the frequency of A is the same as frequency of T, same with C and G pair,
     * and C' and G'.
     *
     * @return frequencies of underlying methylated and unmethylated nucleotides
     */
    @Override
    public double[] getFrequencies() {
        // Used only for root in likelihood function
        // assume equal distribution of A and T, and C and G; and C' and G'
        // stationary distribution depends on the values of kappa, alpha, beta and gamma, not on frequencies of 0, 1 and W in data
        return stationary_distribution();
    }


    /** Calculate stationary distribution using SVD method
     *
     * Calculate stationary distribution pi from rate matrix Q by solving the system of linear equation:
     * Ax = b
     * This solution is the solution to pi Q = 0 with condition sum(pi) = 1.
     *
     * This is calculated by first transposing rate matrix Q and adding row corresponding to the sum condition.
     * Ax = b
     * and then using Singular Value Decomposition to find a pseudo-inverse of the matrix A, so that:
     * A^-1 * b = x
     *
     * The matrix A with the vector b has following form:
     *    A   | b
     * ------------
     *   t(Q) | 0
     * 1 .. 1 | 1
     *
     * @param Q rate matrix
     * @return stationary distribution pi
     */
    private double[] stationary_distribution(double[][] Q) {
        double[][] A = new double[7][6];
        double[] b = new double[]{0, 0, 0, 0, 0, 0, 1};
        // transpose and copy rate matrix Q into matrix A
        for(int i=0; i<6; i++){
            for(int j=0; j<6; j++){
                A[i][j] = Q[j][i];
            }
        }
        // solve
        RealMatrix AA = new Array2DRowRealMatrix(A);
        double[] pi = new SingularValueDecomposition(AA).getSolver().getInverse().operate(b);
        return pi;
    }


    /** Calculate stationary distribution analytically
     *
     * This is an analytical solution to pi*Q = 0 with condition sum(pi) = 1 that was calculated manually.
     *
     * pi(C') and pi(G') = 1 / 2*( 2A + B + 1) = C
     * pi(A) and pi(T) = C * (A + B)
     * pi(C) and pi(G) = C * A
     *
     * where:
     * A = (2 + kappa + beta + gamma) / alpha
     * B = (kappa + gamma + 1) / (kappa + 1)
     *
     *
     * @return stationary distribution pi
     */
    private double[] stationary_distribution(){
        double kappa = kappaPar.getValue();
        double alpha = alphaPar.getValue();
        double beta = betaPar.getValue();
        double gamma = gammaPar.getValue();

        double[] pi = new double[6];
        if(alpha == 0){
            // degenerated case
            pi[0] = pi[1] = pi[2] = pi[3] = 0.25;
        } else {

            double A = (2 + kappa + beta + gamma) / alpha;
            double B = (kappa + gamma + 1) / (kappa + 1);
            double C = 1 / (2 * (2 * A + B + 1) );

            pi[4] = pi[5] = C;
            pi[0] = pi[3] = C * (A + B);
            pi[1] = pi[2] = C * A;
        }

        return pi;
    }


    /** Setup rate matrix for MethylationHKY model.
     *
     * MethylationHKY model is based on HKY model, but with few changes.
     * First of all, model is extended for C' and G' (or MetC and MetC on the opposite strand).
     * Secondly, due to the structure, the model is time non-reversible and thus
     * frequencies of stationary distribution are not used. This significantly limits
     * the rates and number of parameters.
     *
     * The rate matrix has four parameters:
     * kappa -- transversion rate
     * alpha -- methylation rate
     * beta -- demethylation rate
     * gamma -- demethylation into T
     *
     * The rate matrix has following form:
     *
     *       A  C G  T | C' G'
     *   A   -  1 k  1 | 0  0
     *   C   1  - 1  k | a  0
     *   G   k  1 -  1 | 0  a
     *   T   1  k 1  - | 0  0
     *      -----------------
     *   C'  1  b 1 k+g  -  0
     *   G' k+g 1 b  1   0  -
     */
    private double[][] setupUnnormalizedRateMatrix(){
        double[][] unnormalizedRateMatrix = new double[6][6];
        double kappa = kappaPar.getValue(0);
        double alpha = alphaPar.getValue(0);
        double beta = betaPar.getValue(0);
        double gamma = gammaPar.getValue(0);

        // set all states to 0, states should be always initialized as 0, so this shouldn't be required
        for(int i = 0; i < 6; i++){
            for(int j = 0; j < 6; j++){
                unnormalizedRateMatrix[i][j] = 0;
            }
        }

        // set up left side to 1 times frequency
        for(int i = 0; i < 6; i++){
            for(int j = 0; j < 4; j++){
                if(i == j){
                    continue;
                }
                unnormalizedRateMatrix[i][j] = 1;
            }
        }

        // fill in kappa and kappa + gamma
        unnormalizedRateMatrix[0][2] = kappa;
        unnormalizedRateMatrix[2][0] = kappa;
        unnormalizedRateMatrix[1][3] = kappa;
        unnormalizedRateMatrix[3][1] = kappa;
        unnormalizedRateMatrix[5][0] = (kappa + gamma);
        unnormalizedRateMatrix[4][3] = (kappa + gamma);
        // set up beta
        unnormalizedRateMatrix[4][1] = beta;
        unnormalizedRateMatrix[5][2] = beta;
        // set up alpha
        unnormalizedRateMatrix[1][4] = alpha;
        unnormalizedRateMatrix[2][5] = alpha;

        // set up diagonal
        for (int i = 0; i < 6; i++) {
            double sum = 0.0;
            for (int j = 0; j < 6; j++) {
                if (i != j) {
                    sum += unnormalizedRateMatrix[i][j];
                }
            }
            unnormalizedRateMatrix[i][i] = -sum;
        }
        return unnormalizedRateMatrix;
    }


    /** Setup normalized rate matrix for MethylationHKY model.
     *
     * MethylationHKY model is based on HKY model, but with few changes.
     * First of all, model is extended for C' and G' (or MetC and MetC on the opposite strand).
     * Secondly, due to the structure, the model is time non-reversible and thus
     * frequencies of stationary distribution are not used. This significantly limits
     * the rates and number of parameters.
     *
     * The rate matrix has four parameters:
     * kappa -- transversion rate
     * alpha -- methylation rate
     * beta -- demethylation rate
     * gamma -- demethylation into T
     *
     * The rate matrix has following form:
     *
     *       A  C G  T | C' G'
     *   A   -  1 k  1 | 0  0
     *   C   1  - 1  k | a  0
     *   G   k  1 -  1 | 0  a
     *   T   1  k 1  - | 0  0
     *      -----------------
     *   C'  1  b 1 k+g  -  0
     *   G' k+g 1 b  1   0  -
     *
     *   Matrix is then normalized so that expected rate of change is 1.
     */
    @Override
    protected void setupRateMatrix(){
        double[][] unnormalizedRateMatrix = setupUnnormalizedRateMatrix();

        // Normalize rate matrix to one unit per time:
        normalize(unnormalizedRateMatrix);
    }


    /** Normalization of rate matrix for expected rate to be 1
     *
     * This is done by utilizing the stationary distribution pi and calculating
     * the expected number of substitution mu:
     * mu = pi * diag(Q).
     * The normalized rate matrix Q' is then calculated by dividing Q by mu:
     * Q' = Q/mu.
     * Now, the expected rate of Q' should be 1:
     * pi * diag(Q') = 1.
     *
     */
    private void normalize(double[][] unnormalizedRateMatrix){
        double[] pi = stationary_distribution();

        // get the expected rate of non-normalized matrix
        double mu = 0;
        for (int i = 0; i < 6; i++){
            mu -= pi[i] * unnormalizedRateMatrix[i][i];
        }

        // normalize matrix using mu
        for(int i = 0; i < 6; i++){
            for(int j = 0; j < 6; j++){
                rateMatrix[i][j] = unnormalizedRateMatrix[i][j] / mu;
            }
        }
    }

}
