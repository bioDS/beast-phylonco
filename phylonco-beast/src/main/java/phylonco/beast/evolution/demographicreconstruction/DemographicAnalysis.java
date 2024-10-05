package phylonco.beast.evolution.demographicreconstruction;



import lphy.base.evolution.coalescent.populationmodel.ExponentialPopulation;
import lphy.base.evolution.coalescent.populationmodel.GompertzPopulation_f0;
import lphy.base.evolution.coalescent.populationmodel.GompertzPopulation_t50;
import lphy.base.evolution.coalescent.populationmodel.LogisticPopulation;

import java.io.*;
import java.util.*;

public class DemographicAnalysis {

    public static void main(String[] args) {
        // Specify the population model type: "constant", "exponential", "logistic", "gompertz_f0", or "gompertz_t50"
        String modelType = "gompertz_f0"; // Modify as needed

        // Input and output file paths
        String logFilePath = "/Users/xuyuan/workspace/beast-phylonco/gompCoal-40.log";  // Input file path
        String outputCsvPath = "/Users/xuyuan/Desktop/output_results.csv";  // Output file path

        // Define the parameters to read based on the model type
        String[] parameterNames;
        if (modelType.equals("constant")) {
            parameterNames = new String[]{"N0", "tree.height"};
        } else if (modelType.equals("exponential")) {
            parameterNames = new String[]{"N0", "GrowthRate", "tree.height"};
        } else if (modelType.equals("logistic")) {
            parameterNames = new String[]{"t50", "nCarryingCapacity", "b", "tree.height"};
        } else if (modelType.equals("gompertz_f0")) {
            parameterNames = new String[]{"N0", "f0", "b", "tree.height"};
        } else if (modelType.equals("gompertz_t50")) {
            parameterNames = new String[]{"t50", "b", "NInfinity", "tree.height"};
        } else {
            System.err.println("Unknown model type: " + modelType);
            return;
        }

        // Read parameter samples from the log file
        Map<String, List<Double>> parameterSamples = readParameterSamples(logFilePath, parameterNames);

        // Check if essential parameters exist and have matching sample counts
        if (!checkParameters(modelType, parameterSamples)) {
            return; // Exit if parameters are missing or inconsistent
        }

        int sampleCount = parameterSamples.get("tree.height").size();

        // Extract parameter samples
        List<Double> N0_samples = parameterSamples.get("N0");
        List<Double> growthRate_samples = parameterSamples.get("GrowthRate");
        List<Double> t50_samples = parameterSamples.get("t50");
        List<Double> nCarryingCapacity_samples = parameterSamples.get("nCarryingCapacity");
        List<Double> b_samples = parameterSamples.get("b");
        List<Double> f0_samples = parameterSamples.get("f0");
        List<Double> NInfinity_samples = parameterSamples.get("NInfinity");
        List<Double> treeHeight_samples = parameterSamples.get("tree.height");

        // Compute the 95% HPD interval of tree.height
        double[] hpdInterval = calculateHPD(treeHeight_samples, 0.95);
        double hpdUpper = hpdInterval[1];

        // Set the time range
        double minTime = 0.0;
        double maxTime = hpdUpper;

        // Define the time grid
        int binCount = 100;  // Number of time points
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

        // Compute population sizes at each time point for each sample
        for (int sampleIndex = 0; sampleIndex < sampleCount; sampleIndex++) {
            if (modelType.equals("gompertz_f0")) {
                // GompertzPopulation_f0 model
                double N0 = N0_samples.get(sampleIndex);
                double f0 = f0_samples.get(sampleIndex);
                double b = b_samples.get(sampleIndex);

                GompertzPopulation_f0 gompertzPop = new GompertzPopulation_f0(N0, f0, b);

                for (int i = 0; i < timePoints.size(); i++) {
                    double time = timePoints.get(i);
                    double size = gompertzPop.getTheta(time);
                    populationSizesAtTimePoints.get(i).add(size);
                }
            } else if (modelType.equals("gompertz_t50")) {
                // GompertzPopulation_t50 model
                double t50 = t50_samples.get(sampleIndex);
                double b = b_samples.get(sampleIndex);
                double NInfinity = NInfinity_samples.get(sampleIndex);

                GompertzPopulation_t50 gompertzPop = new GompertzPopulation_t50(t50, b, NInfinity);

                for (int i = 0; i < timePoints.size(); i++) {
                    double time = timePoints.get(i);
                    double size = gompertzPop.getTheta(time);
                    populationSizesAtTimePoints.get(i).add(size);
                }
            } else if (modelType.equals("logistic")) {
                // Logistic model
                double t50 = t50_samples.get(sampleIndex);
                double nCarryingCapacity = nCarryingCapacity_samples.get(sampleIndex);
                double b = b_samples.get(sampleIndex);

                LogisticPopulation logPop = new LogisticPopulation(t50, nCarryingCapacity, b);

                for (int i = 0; i < timePoints.size(); i++) {
                    double time = timePoints.get(i);
                    double size = logPop.getTheta(time);
                    populationSizesAtTimePoints.get(i).add(size);
                }
            } else if (modelType.equals("exponential")) {
                // Exponential model
                double N0 = N0_samples.get(sampleIndex);
                double growthRate = growthRate_samples.get(sampleIndex);

                ExponentialPopulation expPop = new ExponentialPopulation(growthRate, N0);

                for (int i = 0; i < timePoints.size(); i++) {
                    double time = timePoints.get(i);
                    double size = expPop.getTheta(time);
                    populationSizesAtTimePoints.get(i).add(size);
                }
            } else if (modelType.equals("constant")) {
                // Constant model
                double N0 = N0_samples.get(sampleIndex);
                for (int i = 0; i < timePoints.size(); i++) {
                    double size = N0; // Population size remains constant
                    populationSizesAtTimePoints.get(i).add(size);
                }
            } else {
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

    // Method to check if essential parameters exist and have matching sample counts
    public static boolean checkParameters(String modelType, Map<String, List<Double>> parameterSamples) {
        List<Double> treeHeight_samples = parameterSamples.get("tree.height");

        if (treeHeight_samples == null || treeHeight_samples.isEmpty()) {
            System.err.println("tree.height parameter samples not found.");
            return false;
        }

        int sampleCount = treeHeight_samples.size();

        if (modelType.equals("constant")) {
            return checkSampleConsistency(parameterSamples, new String[]{"N0"}, sampleCount);
        } else if (modelType.equals("exponential")) {
            return checkSampleConsistency(parameterSamples, new String[]{"N0", "GrowthRate"}, sampleCount);
        } else if (modelType.equals("logistic")) {
            return checkSampleConsistency(parameterSamples, new String[]{"t50", "nCarryingCapacity", "b"}, sampleCount);
        } else if (modelType.equals("gompertz_f0")) {
            return checkSampleConsistency(parameterSamples, new String[]{"N0", "f0", "b"}, sampleCount);
        } else if (modelType.equals("gompertz_t50")) {
            return checkSampleConsistency(parameterSamples, new String[]{"t50", "b", "NInfinity"}, sampleCount);
        } else {
            System.err.println("Unknown model type: " + modelType);
            return false;
        }
    }

    // Helper method to check if specified parameters exist and have matching sample counts
    public static boolean checkSampleConsistency(Map<String, List<Double>> parameterSamples, String[] paramNames, int sampleCount) {
        for (String paramName : paramNames) {
            List<Double> samples = parameterSamples.get(paramName);
            if (samples == null || samples.isEmpty()) {
                System.err.println(paramName + " parameter samples not found.");
                return false;
            }
            if (samples.size() != sampleCount) {
                System.err.println("Sample counts for " + paramName + " and tree.height do not match.");
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
