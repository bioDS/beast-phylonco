

/**
 * DemographicAnalysis
 *
 * This Java program performs demographic analysis based on various population models,
 * including constant, exponential, logistic, Gompertz with f0 or t50 parameters, and expansion models.
 * It reads parameter samples from a log file, computes population sizes over a time grid, and
 * outputs the results to a CSV file.
 *
 * Options:
 * - modelType: Select between "constant", "exponential", "logistic", "gompertz_f0", "gompertz_t50", or "expansion".
 * - treeHeightMethod: Choose how to compute maxTime using Tree.height ("mean", "median", "lowerHPD", "upperHPD", or "fixed").
 * - binCount: Specify the number of time points for the analysis.
 */




package phylonco.beast.evolution.demographicreconstruction;

import lphy.base.evolution.coalescent.PopulationFunction;
import lphy.base.evolution.coalescent.populationmodel.*;

import java.io.*;
import java.util.*;

public class DemographicAnalysis {

    public static void main(String[] args) {
        // Specify the population model type: "constant", "exponential", "logistic", "gompertz_f0", "gompertz_t50", or "expansion"
        String modelType = "expansion"; // Modify as needed

        // Input and output file paths
        String logFilePath = "/Users/xuyuan/workspace/beast-phylonco/hcvskyline copy.log";  // Input file path
        String outputCsvPath = "/Users/xuyuan/Desktop/output_results.csv";  // Output file path

        // **New** Choose method to set maxTime based on Tree.height (mean, median, lower HPD, upper HPD)
        String treeHeightMethod = "lowerHPD"; // Options: "mean", "median", "lowerHPD", "upperHPD", "fixed"

        // **New** Set the maximum binCount (number of time points)
        int binCount = 1000;  // You can modify this value to set the maximum binCount

        // Define the parameters to read based on the model type
        String[] parameterNames = determineParameterNames(modelType);
        if (parameterNames == null) {
            System.err.println("Unknown model type: " + modelType);
            return;
        }

        // Read parameter samples from the log file
        Map<String, List<Double>> parameterSamples = readParameterSamples(logFilePath, parameterNames);

        // Check if essential parameters exist and have matching sample counts
        if (!checkParameters(modelType, parameterSamples)) {
            return; // Exit if parameters are missing or inconsistent
        }

        // Extract tree height samples to set maxTime
        List<Double> treeHeight_samples = parameterSamples.get("Tree.height");
        if (treeHeight_samples == null || treeHeight_samples.isEmpty()) {
            System.err.println("Error: Tree.height samples not found or empty.");
            return;
        }

        // **New** Compute maxTime based on treeHeightMethod (mean, median, lower HPD, upper HPD)
        double maxTime = 0.0;
        switch (treeHeightMethod) {
            case "mean":
                maxTime = calculateMean(treeHeight_samples);
                System.out.println("Max Time based on mean of Tree.height: " + maxTime);
                break;
            case "median":
                maxTime = calculateMedian(treeHeight_samples);
                System.out.println("Max Time based on median of Tree.height: " + maxTime);
                break;
            case "lowerHPD":
                maxTime = calculateHPD(treeHeight_samples, 0.95)[0];
                System.out.println("Max Time based on lower 95% HPD of Tree.height: " + maxTime);
                break;
            case "upperHPD":
                maxTime = calculateHPD(treeHeight_samples, 0.95)[1];
                System.out.println("Max Time based on upper 95% HPD of Tree.height: " + maxTime);
                break;
            case "fixed":
                maxTime = 100.0; // Example of a fixed maxTime
                System.out.println("Max Time set to a fixed value: " + maxTime);
                break;
            default:
                System.err.println("Unknown tree height method.");
                return;
        }

        // Set the time range
        double minTime = 0.0;

        // **New** Define the time grid based on the user-specified binCount
        double deltaTime = (maxTime - minTime) / (binCount - 1);
        List<Double> timePoints = new ArrayList<>();
        for (int i = 0; i < binCount; i++) {
            double time = minTime + i * deltaTime;
            timePoints.add(time);
        }

        // Initialize population size lists
        List<List<Double>> populationSizesAtTimePoints = new ArrayList<>();
        for (int i = 0; i < timePoints.size(); i++) {
            populationSizesAtTimePoints.add(new ArrayList<>());
        }

        // Extract other parameter samples and use them in different models
        for (int sampleIndex = 0; sampleIndex < treeHeight_samples.size(); sampleIndex++) {
            switch (modelType) {
                case "gompertz_f0":
                    double N0_f0 = parameterSamples.get("N0").get(sampleIndex);
                    double f0 = parameterSamples.get("f0").get(sampleIndex);
                    double b_f0 = parameterSamples.get("b").get(sampleIndex);

                    GompertzPopulation_f0 gompertzPop_f0 = new GompertzPopulation_f0(N0_f0, f0, b_f0);
                    calculatePopulationSizes(gompertzPop_f0, timePoints, populationSizesAtTimePoints, sampleIndex);
                    break;

                case "gompertz_t50":
                    double t50 = parameterSamples.get("t50").get(sampleIndex);
                    double b_t50 = parameterSamples.get("b").get(sampleIndex);
                    double NInfinity = parameterSamples.get("NInfinity").get(sampleIndex);

                    GompertzPopulation_t50 gompertzPop_t50 = new GompertzPopulation_t50(t50, b_t50, NInfinity);
                    calculatePopulationSizes(gompertzPop_t50, timePoints, populationSizesAtTimePoints, sampleIndex);
                    break;

                case "logistic":
                    double t50_log = parameterSamples.get("t50").get(sampleIndex);
                    double nCarryingCapacity = parameterSamples.get("nCarryingCapacity").get(sampleIndex);
                    double b_log = parameterSamples.get("b").get(sampleIndex);

                    LogisticPopulation logPop = new LogisticPopulation(t50_log, nCarryingCapacity, b_log);
                    calculatePopulationSizes(logPop, timePoints, populationSizesAtTimePoints, sampleIndex);
                    break;

                case "exponential":
                    double N0_exp = parameterSamples.get("N0").get(sampleIndex);
                    double growthRate = parameterSamples.get("GrowthRate").get(sampleIndex);

                    ExponentialPopulation expPop = new ExponentialPopulation(growthRate, N0_exp);
                    calculatePopulationSizes(expPop, timePoints, populationSizesAtTimePoints, sampleIndex);
                    break;

                case "constant":
                    double N0_const = parameterSamples.get("N0").get(sampleIndex);
                    for (int i = 0; i < timePoints.size(); i++) {
                        double size = N0_const; // Population size remains constant
                        populationSizesAtTimePoints.get(i).add(size);
                    }
                    break;

                case "expansion":
                    double N0_expansion = parameterSamples.get("N0").get(sampleIndex);
                    double tau = parameterSamples.get("tau").get(sampleIndex);
                    double r = parameterSamples.get("r").get(sampleIndex);
                    double NC = parameterSamples.get("NC").get(sampleIndex);

                    ExpansionPopulation expanPop = new ExpansionPopulation(N0_expansion, tau, r, NC);
                    calculatePopulationSizes(expanPop, timePoints, populationSizesAtTimePoints, sampleIndex);
                    break;

                default:
                    System.err.println("Unknown model type: " + modelType);
                    return;
            }
        }

        // Compute statistics
        List<Double> means = new ArrayList<>();
        List<Double> medians = new ArrayList<>();
        List<Double> lower95 = new ArrayList<>();
        List<Double> upper95 = new ArrayList<>();

        for (int i = 0; i < timePoints.size(); i++) {
            List<Double> sizes = populationSizesAtTimePoints.get(i);

            double mean = calculateMean(sizes);
            double median = calculateMedian(sizes);
            double lower = calculatePercentile(sizes, 2.5);
            double upper = calculatePercentile(sizes, 97.5);

            means.add(mean);
            medians.add(median);
            lower95.add(lower);
            upper95.add(upper);
        }

        // Write results to CSV file
        writeResultsToCSV(outputCsvPath, timePoints, lower95, upper95, means, medians);

        System.out.println("Analysis complete. Results saved to " + outputCsvPath);
    }

    // Method to determine parameter names based on model type
    public static String[] determineParameterNames(String modelType) {
        switch (modelType) {
            case "constant":
                return new String[]{"N0", "Tree.height"};
            case "exponential":
                return new String[]{"N0", "GrowthRate", "Tree.height"};
            case "logistic":
                return new String[]{"t50", "nCarryingCapacity", "b", "Tree.height"};
            case "gompertz_f0":
                return new String[]{"N0", "f0", "b", "Tree.height"};
            case "gompertz_t50":
                return new String[]{"t50", "b", "NInfinity", "Tree.height"};
            case "expansion":
                return new String[]{"N0", "tau", "r", "NC", "Tree.height"};
            default:
                return null;
        }
    }

    // Method to calculate population sizes
    public static void calculatePopulationSizes(PopulationFunction populationFunction, List<Double> timePoints, List<List<Double>> populationSizesAtTimePoints, int sampleIndex) {
        for (int i = 0; timePoints.size() > i; i++) {
            double time = timePoints.get(i);
            double size = populationFunction.getTheta(time);
            populationSizesAtTimePoints.get(i).add(size);
        }
    }

    // Method to check if essential parameters exist and have matching sample counts
    public static boolean checkParameters(String modelType, Map<String, List<Double>> parameterSamples) {
        List<Double> treeHeight_samples = parameterSamples.get("Tree.height");

        if (treeHeight_samples == null || treeHeight_samples.isEmpty()) {
            System.err.println("Tree.height parameter samples not found.");
            return false;
        }

        int sampleCount = treeHeight_samples.size();

        String[] requiredParams = determineParameterNames(modelType);
        if (requiredParams == null) {
            System.err.println("Unknown model type: " + modelType);
            return false;
        }

        for (String param : requiredParams) {
            List<Double> samples = parameterSamples.get(param);
            if (samples == null || samples.isEmpty()) {
                System.err.println(param + " parameter samples not found.");
                return false;
            }
            if (samples.size() != sampleCount) {
                System.err.println("Sample counts for " + param + " and Tree.height do not match.");
                return false;
            }
        }

        return true;
    }

    // Method to read parameter samples from the log file
    public static Map<String, List<Double>> readParameterSamples(String filePath, String[] parameterNames) {
        Map<String, List<Double>> samples = new HashMap<>();
        for (String param : parameterNames) {
            samples.put(param, new ArrayList<>());
        }
        try (BufferedReader br = new BufferedReader(new FileReader(filePath))) {
            String line;
            String[] headers = null;

            // Detect delimiter (tab or whitespace)
            String delimiter = "\\s+"; // Default to whitespace

            while ((line = br.readLine()) != null) {
                line = line.trim();
                // Skip empty lines or comments
                if (line.isEmpty() || line.startsWith("#")) {
                    continue;
                }
                // Handle header line
                if (headers == null) {
                    if (line.contains("\t")) {
                        delimiter = "\t"; // Use tab as delimiter
                    }
                    headers = line.split(delimiter);
                    continue;
                }
                // Read data lines
                String[] tokens = line.split(delimiter);
                if (tokens.length != headers.length) {
                    continue; // Skip malformed lines
                }
                Map<String, String> dataMap = new HashMap<>();
                for (int i = 0; i < headers.length; i++) {
                    dataMap.put(headers[i], tokens[i]);
                }
                // Extract required parameters
                for (String param : parameterNames) {
                    String valueStr = dataMap.get(param);
                    if (valueStr != null) {
                        double value = Double.parseDouble(valueStr);
                        samples.get(param).add(value);
                    }
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        return samples;
    }

    // Compute mean
    public static double calculateMean(List<Double> values) {
        double sum = 0.0;
        for (double val : values) {
            sum += val;
        }
        return sum / values.size();
    }

    // Compute median
    public static double calculateMedian(List<Double> values) {
        List<Double> sorted = new ArrayList<>(values);
        Collections.sort(sorted);
        int n = sorted.size();
        if (n % 2 == 0) {
            return (sorted.get(n / 2 - 1) + sorted.get(n / 2)) / 2.0;
        } else {
            return sorted.get(n / 2);
        }
    }

    // Compute a specific percentile
    public static double calculatePercentile(List<Double> values, double percentile) {
        List<Double> sorted = new ArrayList<>(values);
        Collections.sort(sorted);
        int index = (int) Math.ceil(percentile / 100.0 * sorted.size()) - 1;
        index = Math.max(0, index);
        index = Math.min(index, sorted.size() - 1);
        return sorted.get(index);
    }

    // Compute HPD interval
    public static double[] calculateHPD(List<Double> values, double credMass) {
        List<Double> sorted = new ArrayList<>(values);
        Collections.sort(sorted);
        int n = sorted.size();
        int intervalIdxInc = (int) Math.floor(credMass * n);
        double minWidth = Double.MAX_VALUE;
        int hpdStartIdx = 0;

        for (int i = 0; i <= (n - intervalIdxInc); i++) {
            double width = sorted.get(i + intervalIdxInc - 1) - sorted.get(i);
            if (width < minWidth) {
                minWidth = width;
                hpdStartIdx = i;
            }
        }
        double hpdLower = sorted.get(hpdStartIdx);
        double hpdUpper = sorted.get(hpdStartIdx + intervalIdxInc - 1);
        return new double[]{hpdLower, hpdUpper};
    }

    // Write results to CSV file
    public static void writeResultsToCSV(String filePath, List<Double> timePoints, List<Double> lower95, List<Double> upper95, List<Double> means, List<Double> medians) {
        try (PrintWriter pw = new PrintWriter(new FileWriter(filePath))) {
            pw.println("Time,Lower95,Upper95,Mean,Median");
            for (int i = 0; i < timePoints.size(); i++) {
                pw.printf("%f,%f,%f,%f,%f%n", timePoints.get(i), lower95.get(i), upper95.get(i), means.get(i), medians.get(i));
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
