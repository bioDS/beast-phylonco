package phylonco.lphy.evolution.copynumbermodel;

import org.junit.Test;

import lphy.core.model.Value;

import static org.junit.Assert.assertTrue;

public class NestedBDTest {
    /**
     *
     */

//For Birth-Death process
    public void simulateCopiesOnBranchBinTest(boolean debug, boolean unconditional, int trials, int numSim) {
        Value<Double> lambda = new Value<>("lambda", 0.7);
        Value<Double> mu = new Value<>("mu", 0.3);
        Value<Integer> nBins = new Value<>("nBins", 5);

        CopyNumberBD model = new CopyNumberBD(lambda, mu, nBins);
        double startingCopies = 2;
        double time = 3;

        // Calculate Theoretical Unconditional (including extinctions) mean
        double netRate = lambda.value() - mu.value();
        double expTerm = Math.exp(netRate * time);
        double meanExpUncond = startingCopies * expTerm;
        // Calculate extinction probability
        double q;
        if (Math.abs(netRate) > 1e-10) { // λ ≠ μ
            double qTerm = (mu.value() * (expTerm - 1)) / (lambda.value() * expTerm - mu.value());
            q = Math.pow(qTerm, startingCopies);  // Apply power for initial population size
        } else { // λ = μ
            double qTerm = mu.value() * time / (1 + mu.value() * time);
            q = Math.pow(qTerm, startingCopies);
        }
        double survivalProb = 1 - q;
        // Calculate Theoretical Conditional (excluding extinctions) expected mean
        double meanExpCond = meanExpUncond / survivalProb;

        // Simulations
        for (int trial = 0; trial < trials; trial++) {
            int[] results = new int[numSim];
            int sum = 0;
            int extinctCount = 0;
            // Run multiple simulations for each trial
            for (int i = 0; i < numSim; i++) {
                results[i] = model.simulateCopiesOnBranchBin(0, time, (int) startingCopies);
                sum += results[i];
                if (results[i] == 0) {
                    extinctCount++;
                }
            }

            // Calculate unconditional mean
            double meanObsUncond = sum / (double) numSim;
            // Calculate unconditional
            double sumSquaredDevsUncond = 0;
            for (int i = 0; i < numSim; i++) {
                double deviation = results[i] - meanObsUncond;
                sumSquaredDevsUncond += deviation * deviation;
            }
            double varObsUncond = sumSquaredDevsUncond / (numSim - 1);
            double stdErrorUncond = Math.sqrt(varObsUncond / numSim);
            // Calculate the Unconditional expected range (95% CI)
            double lowerBoundUncond = meanObsUncond - 1.96 * stdErrorUncond;
            double upperBoundUncond = meanObsUncond + 1.96 * stdErrorUncond;
            // Check if Unconditional expected mean falls within the range
            boolean withinRangeUncond = (meanExpUncond >= lowerBoundUncond && meanExpUncond <= upperBoundUncond);

            // Calculate conditional mean
            int trialsWithoutZero = numSim - extinctCount;
            if (extinctCount < numSim) { // Prevent division by zero
                double meanObsCond = (double) sum / trialsWithoutZero;
                // Calculate conditional variance and stdError
                double sumSquaredDevsCond = 0;
                for (int i = 0; i < numSim; i++) {
                    if (results[i] > 0) {
                        double deviation = results[i] - meanObsCond;
                        sumSquaredDevsCond += deviation * deviation;
                    }
                }
                double varObsCond = sumSquaredDevsCond / (trialsWithoutZero - 1);
                double stdErrorCond = Math.sqrt(varObsCond / trialsWithoutZero);
                // Calculate the Conditional expected range using standard deviation (95% CI)
                double lowerBoundCond = meanObsCond - 1.96 * stdErrorCond;
                double upperBoundCond = meanObsCond + 1.96 * stdErrorCond;
                // Check if Conditional expected mean falls within the range
                boolean withinRangeCond = (meanExpCond >= lowerBoundCond && meanExpCond <= upperBoundCond);


                // csv format
                if (debug) {
                    if (unconditional) {
                        // Print unconditional results only
                        System.out.print(meanExpUncond + ",");
                        System.out.print(meanObsUncond + ",");
                        System.out.print(lowerBoundUncond + "," + upperBoundUncond + ",");
                        System.out.println(withinRangeUncond);
                    } else if (extinctCount < numSim) { // Only if there are non-extinct cases
                        // Print conditional results only
                        System.out.print(meanExpCond + ",");
                        System.out.print(meanObsCond + ",");
                        System.out.print(lowerBoundCond + "," + upperBoundCond + ",");
                        System.out.println(withinRangeCond);
                    }
                } else {
                    // Test assertions for both conditions
                    assertTrue(withinRangeCond);
                    assertTrue(withinRangeUncond);
                }
            }
        }
    }

    // Debugging
    @Test
    public void testBirthAndDeathMulti() {
        String header = "meanExp,meanObs,lowerBound,upperBound,withinRange";
        // Print conditional header once
        System.out.println("====== CONDITIONAL (only non-extinct) ======");
        System.out.println(header);
        // Run conditional analysis
        simulateCopiesOnBranchBinTest(true, false,1000,1000);
        // Print unconditional header once
        System.out.println("====== UNCONDITIONAL (including extinction) ======");
        System.out.println(header);
        // Run unconditional analysis
        simulateCopiesOnBranchBinTest(true, true,1000,1000);
    }

    // Test
    @Test
    public void testBirthAndDeathSingle() {
        // Run a single trial with assertions (no output)
        simulateCopiesOnBranchBinTest(false, false,1,1000);
    }


// For Birth-only process
    public void testBirthOnly(boolean debug, boolean unconditional, int trials, int numSim) {
        Value<Double> lambda = new Value<>("lambda", 0.7);
        Value<Double> mu = new Value<>("mu", 0.0);
        Value<Integer> nBins = new Value<>("nBins", 5);

        CopyNumberBD model = new CopyNumberBD(lambda, mu, nBins);
        double startingCopies = 2;
        double time = 3;

        // Calculate theoretical expected mean
        // Conditional and unconditional means are the same in a pure birth process
        // since extinction is impossible
        double expTerm = Math.exp(lambda.value() * time);
        double meanExp = startingCopies * expTerm;

        // Simulations
        for (int trial = 0; trial < trials; trial++) {
            int[] results = new int[numSim];
            int sum = 0;
            // Run multiple simulations for each trial
            for (int i = 0; i < numSim; i++) {
                results[i] = model.simulateCopiesOnBranchBin(0, time, (int) startingCopies);
                sum += results[i];
                if (results[i] == 0) {
                }
            }

            // Calculate mean
            double meanObs = sum / (double) numSim;
            // Calculate variance and stdError
            double sumSquaredDevs = 0;
            for (int i = 0; i < numSim; i++) {
                double deviation = results[i] - meanObs;
                sumSquaredDevs += deviation * deviation;
            }
            double varObs = sumSquaredDevs / (numSim - 1);
            double stdError = Math.sqrt(varObs / numSim);
            // Calculate the expected range (95% CI)
            double lowerBound = meanObs - 1.96 * stdError;
            double upperBound = meanObs + 1.96 * stdError;
            // Check if expected mean falls within the range
            boolean withinRange = (meanExp >= lowerBound && meanExp <= upperBound);

            // csv format
            if (debug) {
                System.out.print(meanExp + ",");
                System.out.print(meanObs + ",");
                System.out.print(lowerBound + "," + upperBound + ",");
                System.out.println(withinRange);
            } else {
                // Test assertion
                assertTrue(withinRange);
            }
        }
    }
    // Debugging
    @Test
    public void testBirthOnlyMulti() {
        String header = "meanExp,meanObs,lowerBound,upperBound,withinRange";
        // Print header once
        System.out.println("====== PURE BIRTH PROCESS ======");
        System.out.println(header);
        // Run analysis
        testBirthOnly(true, false, 1000, 1000);
    }
    // Test
    @Test
    public void testBirthOnlySingle() {
        // Run a single trial with assertions (no output)
        testBirthOnly(false, false, 1, 1000);
    }
}