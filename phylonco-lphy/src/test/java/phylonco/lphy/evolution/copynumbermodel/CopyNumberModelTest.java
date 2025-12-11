package phylonco.lphy.evolution.copynumbermodel;

import org.junit.Test;
import lphy.core.model.Value;
import static org.junit.Assert.assertTrue;

/**
 * Tests for the copy number birth-death model in {@link CopyNumberBD}.
 *
 * <p>Validates that simulateCopiesOnBranchBin produces statistically correct results
 * by comparing observed means from simulations against theoretical means (within
 * confidence intervals).</p>
 *
 * <h3>Test Scenarios:</h3>
 * <ul>
 * <li><b>Birth-death processes</b> (λ > 0, μ > 0): Can lead to extinction (state 0)</li>
 * <li><b>Pure birth processes</b> (λ > 0, μ = 0): No extinction possible</li>
 * <li><b>Both parameterization methods</b>: separate birth (lambda)/death (mu) rates OR single bdRate (assume λ = μ)</li>
 * </ul>
 *
 * <h3>Conditional vs. Unconditional Expectations:</h3>
 * <p>For birth-death processes where extinction is possible, we test two types of expectations:</p>
 * <ul>
 * <li><b>Unconditional</b>: Average across ALL simulations (including extinctions where state = 0)</li>
 * <li><b>Conditional</b>: Average only among surviving lineages (excluding extinctions)</li>
 * </ul>
 * <p><i>Note:</i> For pure birth processes (μ = 0), extinction is impossible, so
 * conditional and unconditional expectations are identical.</p>
 */
public class CopyNumberModelTest {

    // ==================== Birth-Death Process Tests ====================

    /**
     * Single tests
     */
    @Test
    public void testBirthDeathWithLambdaMu_DifferentRates() {
        // Test with lambda != mu (only works with lambda/mu constructor)
        testBirthDeathProcess(false, false, 1, 1000, 0.7, 0.3, 2);
    }

    @Test
    public void testBirthDeathWithBdRate_EqualRates() {
        // Test with lambda = mu (using bdRate constructor)
        testBirthDeathProcessWithBdRate(false, false, 1, 1000, 0.5, 2);
    }

    /**
     * Multi-tests with Different birth-death Rates (Conditional and Unconditional)
     */
    @Test
    public void testBirthDeathConditionalMulti_LambdaMu() {
        System.out.println("====== BIRTH-DEATH CONDITIONAL (lambda/mu, lambda != mu) ======");
        System.out.println("meanExp,meanObs,lowerBound,upperBound,withinRange");
        testBirthDeathProcess(true, false, 1000, 1000, 0.7, 0.3, 2);
    }

    @Test
    public void testBirthDeathUnconditionalMulti_LambdaMu() {
        System.out.println("====== BIRTH-DEATH UNCONDITIONAL (lambda/mu, lambda != mu) ======");
        System.out.println("meanExp,meanObs,lowerBound,upperBound,withinRange");
        testBirthDeathProcess(true, true, 1000, 1000, 0.7, 0.3, 2);
    }

    /**
     * Multi-tests with Equal birth-death Rates (Conditional and Unconditional)
     */
    @Test
    public void testBirthDeathConditionalMulti_BdRate() {
        System.out.println("====== BIRTH-DEATH CONDITIONAL (bdRate, lambda = mu) ======");
        System.out.println("meanExp,meanObs,lowerBound,upperBound,withinRange");
        testBirthDeathProcessWithBdRate(true, false, 1000, 1000, 0.8, 2);
    }

    @Test
    public void testBirthDeathUnconditionalMulti_BdRate() {
        System.out.println("====== BIRTH-DEATH UNCONDITIONAL (bdRate, lambda = mu) ======");
        System.out.println("meanExp,meanObs,lowerBound,upperBound,withinRange");
        testBirthDeathProcessWithBdRate(true, true, 1000, 1000, 0.8, 2);
    }

    // ==================== Pure Birth Process Tests ====================
    //For pure birth, conditional = unconditional

    /**
     * Single tests
     */
    @Test
    public void testPureBirthWithLambdaMu() {
        // Pure birth with lambda/mu constructor (mu = 0)
        testPureBirthProcess(false, 1, 1000, 0.7, 2);
    }
    /**
     * Multi-tests
     */
    @Test
    public void testPureBirthMulti_LambdaMu() {
        System.out.println("====== PURE BIRTH PROCESS (lambda/mu parameterization) ======");
        System.out.println("meanExp,meanObs,lowerBound,upperBound,withinRange");
        testPureBirthProcess(true, 1000, 1000, 0.8, 2);
    }

    // ==================== Core Testing Logic ====================

    /**
     * Test birth-death process using lambda/mu constructor.
     */
    private void testBirthDeathProcess(boolean debug, boolean unconditional,
                                       int trials, int numSim,
                                       double lambda, double mu, int rootState) {
        testBirthDeathProcessInternal(debug, unconditional, trials, numSim,
                false, lambda, mu, rootState);
    }

    /**
     * Test birth-death process using bdRate constructor (lambda = mu).
     */
    private void testBirthDeathProcessWithBdRate(boolean debug, boolean unconditional,
                                                 int trials, int numSim,
                                                 double bdRate, int rootState) {
        testBirthDeathProcessInternal(debug, unconditional, trials, numSim,
                true, bdRate, bdRate, rootState);
    }

    /**
     * Internal method for birth-death process testing.
     *
     * @param debug         Print CSV output
     * @param unconditional Test unconditional (true) or conditional (false) expectations
     * @param trials        Number of independent trials
     * @param numSim        Number of simulations per trial
     * @param useBdRate     Use bdRate constructor (true) or lambda/mu constructor (false)
     * @param lambda     Birth rate
     * @param mu         Death rate
     * @param rootState  Initial copy number state at root
     */
    private void testBirthDeathProcessInternal(boolean debug, boolean unconditional,
                                               int trials, int numSim, boolean useBdRate,
                                               double lambda, double mu, int rootState) {

        // Create model using appropriate constructor
        CopyNumberBD model = createModel(lambda, mu, rootState, useBdRate);

        // Simulation parameters
        double time = 3.0;
        int startingCopies = rootState;

        // Calculate theoretical expectations
        double netRate = lambda - mu;
        double expTerm = Math.exp(netRate * time);
        double meanExpUncond = startingCopies * expTerm;

        // Calculate extinction probability
        double q = calculateExtinctionProbability(lambda, mu, time, startingCopies, netRate, expTerm);
        double survivalProb = 1 - q;
        double meanExpCond = meanExpUncond / survivalProb;

        // Run trials
        for (int trial = 0; trial < trials; trial++) {
            SimulationResults results = runSimulations(model, time, startingCopies, numSim);

            if (unconditional) {
                checkExpectation(debug, meanExpUncond, results.meanUncond,
                        results.stdErrorUncond, "unconditional");
            } else {
                if (results.nonExtinctCount > 0) {
                    checkExpectation(debug, meanExpCond, results.meanCond,
                            results.stdErrorCond, "conditional");
                }
            }
        }
    }

    /**
     * Test pure birth process using lambda/mu constructor (mu = 0, extinction impossible)).
     */
    private void testPureBirthProcess(boolean debug, int trials, int numSim,
                                      double lambda, int rootState) {
        // Model parameters
        double mu = 0.0;

        // Create model using appropriate constructor
        CopyNumberBD model = createModel(lambda, mu, rootState, false);

        // Simulation parameters
        double time = 3.0;
        int startingCopies = rootState;

        // For pure birth, conditional = unconditional
        double expTerm = Math.exp(lambda * time);
        double meanExp = startingCopies * expTerm;

        // Run trials
        for (int trial = 0; trial < trials; trial++) {
            SimulationResults results = runSimulations(model, time, startingCopies, numSim);
            checkExpectation(debug, meanExp, results.meanUncond,
                    results.stdErrorUncond, "pure birth");
        }
    }

    // ==================== Helper Methods ====================

    /**
     * Create CopyNumberBD model using either constructor
     */
    private CopyNumberBD createModel(double lambda, double mu, int rootState, boolean useBdRate) {
        Value<Integer> rootStateValue = new Value<>("rootState", rootState);
        Value<Integer> nstateValue =null; //Use default (15) by passing null

        if (useBdRate) {
            // Use bdRate constructor (only valid when lambda == mu)
            if (Math.abs(lambda - mu) < 1e-10) {
                Value<Double> bdRate = new Value<>("bdRate", lambda);
                return new CopyNumberBD(bdRate, rootStateValue,nstateValue);
            } else {
                throw new IllegalArgumentException(
                        "Cannot use bdRate constructor when lambda != mu. Use lambda/mu constructor instead.");
            }
        } else {
            // Use lambda/mu constructor
            Value<Double> lambdaValue = new Value<>("lambda", lambda);
            Value<Double> muValue = new Value<>("mu", mu);
            return new CopyNumberBD(lambdaValue, muValue, rootStateValue,nstateValue);
        }
    }

    /**
     * Calculate extinction probability for birth-death process
     */
    private double calculateExtinctionProbability(double lambda, double mu, double time,
                                                  int startingCopies, double netRate, double expTerm) {
        double q;
        if (Math.abs(netRate) > 1e-10) { // λ ≠ μ
            double qTerm = (mu * (expTerm - 1)) / (lambda * expTerm - mu);
            q = Math.pow(qTerm, startingCopies);
        } else { // λ = μ
            double qTerm = mu * time / (1 + mu * time);
            q = Math.pow(qTerm, startingCopies);
        }
        return q;
    }

    /**
     * Run simulations and collect statistics
     */
    private SimulationResults runSimulations(CopyNumberBD model, double time,
                                             int startingCopies, int numSim) {
        int[] results = new int[numSim];
        int sum = 0;
        int extinctCount = 0;

        // Run simulations
        for (int i = 0; i < numSim; i++) {
            results[i] = model.simulateCopiesOnBranchBin(0, time, startingCopies);
            sum += results[i];
            if (results[i] == 0) {
                extinctCount++;
            }
        }

        // Calculate unconditional statistics
        double meanUncond = sum / (double) numSim;
        double stdErrorUncond = calculateStandardError(results, meanUncond, numSim);

        // Calculate conditional statistics (excluding extinctions)
        int nonExtinctCount = numSim - extinctCount;
        double meanCond = 0.0;
        double stdErrorCond = 0.0;

        if (nonExtinctCount > 0) {
            meanCond = sum / (double) nonExtinctCount;
            stdErrorCond = calculateStandardErrorConditional(results, meanCond, nonExtinctCount);
        }

        return new SimulationResults(meanUncond, stdErrorUncond, meanCond,
                stdErrorCond, nonExtinctCount);
    }

    /**
     * Calculate standard error for all simulations
     */
    private double calculateStandardError(int[] results, double mean, int n) {
        double sumSquaredDevs = 0;
        for (int result : results) {
            double deviation = result - mean;
            sumSquaredDevs += deviation * deviation;
        }
        double variance = sumSquaredDevs / (n - 1);
        return Math.sqrt(variance / n);
    }

    /**
     * Calculate standard error excluding extinctions
     */
    private double calculateStandardErrorConditional(int[] results, double mean, int nonExtinctCount) {
        double sumSquaredDevs = 0;
        for (int result : results) {
            if (result > 0) {
                double deviation = result - mean;
                sumSquaredDevs += deviation * deviation;
            }
        }
        double variance = sumSquaredDevs / (nonExtinctCount - 1);
        return Math.sqrt(variance / nonExtinctCount);
    }

    /**
     * Check if theoretical expectation falls within 95% CI of observed mean
     */
    private void checkExpectation(boolean debug, double meanExp, double meanObs,
                                  double stdError, String testType) {
        double lowerBound = meanObs - 1.96 * stdError;
        double upperBound = meanObs + 1.96 * stdError;
        boolean withinRange = (meanExp >= lowerBound && meanExp <= upperBound);

        if (debug) {
            // CSV format for analysis
            System.out.printf("%.6f,%.6f,%.6f,%.6f,%b%n",
                    meanExp, meanObs, lowerBound, upperBound, withinRange);
        } else {
            // Assertion for unit testing
            assertTrue(String.format(
                            "%s test failed: expected %.4f not in CI [%.4f, %.4f]",
                            testType, meanExp, lowerBound, upperBound),
                    withinRange);
        }
    }

    // ==================== Data Classes ====================

    /**
     * Container for simulation results
     */
    private static class SimulationResults {
        final double meanUncond;
        final double stdErrorUncond;
        final double meanCond;
        final double stdErrorCond;
        final int nonExtinctCount;

        SimulationResults(double meanUncond, double stdErrorUncond,
                          double meanCond, double stdErrorCond, int nonExtinctCount) {
            this.meanUncond = meanUncond;
            this.stdErrorUncond = stdErrorUncond;
            this.meanCond = meanCond;
            this.stdErrorCond = stdErrorCond;
            this.nonExtinctCount = nonExtinctCount;
        }
    }
}