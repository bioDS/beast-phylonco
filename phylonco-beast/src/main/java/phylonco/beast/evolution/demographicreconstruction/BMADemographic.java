package phylonco.beast.evolution.demographicreconstruction;

import lphy.base.evolution.coalescent.PopulationFunction;
import lphy.base.evolution.coalescent.populationmodel.*;

import java.io.*;
import java.util.*;

public class BMADemographic {

    public static void main(String[] args) {
        String modelType = "SVS"; // Using the SVS model which dynamically selects the population function

        // Input and output file paths
        String logFilePath = "/Users/xuyuan/workspace/py/SVS/sim1/beast/SVS_Coal-88.log";  // Input file path
        String outputCsvPath = "/Users/xuyuan/Desktop/output_results.csv";  // Output file path

        System.out.println("Starting SVSAnalysis...");
        System.out.println("Reading file: " + logFilePath);

        // Define the parameters from the log file
        String[] parameterNames = {"I", "N0_con", "b_exp", "N0_exp", "b_log", "Ninf_log", "t50_log",
                "NInf_gomp", "b_gom", "t50_gom", "N0_expansion", "NC_expansion",
                "r_expansion", "tau_expansion", "tree.height"};

        try {
            // Read the parameter samples
            Map<String, List<Double>> parameterSamples = readParameterSamples(logFilePath, parameterNames);
            System.out.println("Parameters loaded successfully.");

            // Extract relevant parameters
            List<Double> I_samples = parameterSamples.get("I");
            List<Double> N0_con_samples = parameterSamples.get("N0_con");
            List<Double> b_exp_samples = parameterSamples.get("b_exp");
            List<Double> N0_exp_samples = parameterSamples.get("N0_exp");
            List<Double> b_log_samples = parameterSamples.get("b_log");
            List<Double> Ninf_log_samples = parameterSamples.get("Ninf_log");
            List<Double> t50_log_samples = parameterSamples.get("t50_log");
            List<Double> NInf_gomp_samples = parameterSamples.get("NInf_gomp");
            List<Double> b_gom_samples = parameterSamples.get("b_gom");
            List<Double> t50_gom_samples = parameterSamples.get("t50_gom");
            List<Double> N0_expansion_samples = parameterSamples.get("N0_expansion");
            List<Double> NC_expansion_samples = parameterSamples.get("NC_expansion");
            List<Double> r_expansion_samples = parameterSamples.get("r_expansion");
            List<Double> tau_expansion_samples = parameterSamples.get("tau_expansion");

            // Time points for population size evaluation
            List<Double> timePoints = generateTimePoints(100, 0.0, parameterSamples.get("tree.height").get(0));

            // Population sizes at each time point (for output)
            List<List<Double>> populationSizesAtTimePoints = new ArrayList<>();
            for (int i = 0; i < timePoints.size(); i++) {
                populationSizesAtTimePoints.add(new ArrayList<>());
            }

            // Loop over samples (using indicator I to select the population model)
            for (int sampleIndex = 0; sampleIndex < I_samples.size(); sampleIndex++) {
                double I = I_samples.get(sampleIndex);
                //System.out.println("Sample Index: " + sampleIndex + ", I value: " + I);

                // Dynamically select the population model based on the value of I
                PopulationFunction selectedModel;
                if (I == 0) {
                    selectedModel = new ConstantPopulation(N0_con_samples.get(sampleIndex));
                } else if (I == 1) {
                    selectedModel = new ExponentialPopulation(b_exp_samples.get(sampleIndex), N0_exp_samples.get(sampleIndex));
                } else if (I == 2) {
                    selectedModel = new LogisticPopulation(t50_log_samples.get(sampleIndex), Ninf_log_samples.get(sampleIndex), b_log_samples.get(sampleIndex));
                } else if (I == 3) {
                    selectedModel = new GompertzPopulation_f0(NInf_gomp_samples.get(sampleIndex), b_gom_samples.get(sampleIndex), t50_gom_samples.get(sampleIndex));
                } else if (I == 4) {
                    selectedModel = new ExpansionPopulation(N0_expansion_samples.get(sampleIndex), tau_expansion_samples.get(sampleIndex), r_expansion_samples.get(sampleIndex), NC_expansion_samples.get(sampleIndex));
                } else {
                    throw new IllegalArgumentException("Invalid model indicator I: " + I);
                }

                // Now we pass the selected model to SVSPopulation
                SVSPopulation svsPop = new SVSPopulation(selectedModel);

                // Calculate population sizes at each time point
                for (int i = 0; i < timePoints.size(); i++) {
                    double time = timePoints.get(i);
                    double populationSize = svsPop.getTheta(time);
                    populationSizesAtTimePoints.get(i).add(populationSize);
                }
            }

            // Compute statistics (mean, median, 95% HPD interval)
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

            // Write results to CSV (including time points, mean, median, and 95% HPD intervals)
            System.out.println("Writing results to CSV: " + outputCsvPath);
            writeResultsToCSV(outputCsvPath, timePoints, lower95, upper95, means, medians);

        } catch (Exception e) {
            System.err.println("An error occurred during the analysis.");
            e.printStackTrace();
        }

        System.out.println("Analysis complete.");
    }

    // Generate time points for population size evaluation
    public static List<Double> generateTimePoints(int numPoints, double minTime, double maxTime) {
        List<Double> timePoints = new ArrayList<>();
        double deltaTime = (maxTime - minTime) / (numPoints - 1);
        for (int i = 0; i < numPoints; i++) {
            double time = minTime + i * deltaTime;
            timePoints.add(time);
        }
        return timePoints;
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
            System.err.println("Error reading the file: " + filePath);
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

    // Write results to CSV file (including mean, median, and 95% HPD intervals)
    public static void writeResultsToCSV(String filePath, List<Double> timePoints, List<Double> lower95, List<Double> upper95, List<Double> means, List<Double> medians) {
        try (PrintWriter pw = new PrintWriter(new FileWriter(filePath))) {
            pw.println("Time,Lower95,Upper95,Mean,Median");
            for (int i = 0; i < timePoints.size(); i++) {
                pw.printf("%f,%f,%f,%f,%f%n", timePoints.get(i), lower95.get(i), upper95.get(i), means.get(i), medians.get(i));
            }
        } catch (IOException e) {
            System.err.println("Error writing to CSV file: " + filePath);
            e.printStackTrace();
        }
    }
}
