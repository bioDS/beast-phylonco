package phylonco.lphy.evolution.copynumbermodel;

import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Compresses CNV data by merging consecutive rows with identical copy number patterns
 * Example input file: {@code examples/data/copyNumberSim.txt}
 * Example output file: {@code examples/data/copyNumberSim_compressed.txt}
 */

public class CopyNumberCompressor {

    /**
     * Represents one row of CNV data (one genomic bin)
     */
    public static class CNVRow {
        public String chromosome;
        public long start;
        public long end;
        public int[] copyNumbers; // Copy numbers for all cells

        public CNVRow(String chromosome, long start, long end, int[] copyNumbers) {
            this.chromosome = chromosome;
            this.start = start;
            this.end = end;
            this.copyNumbers = copyNumbers.clone();
        }

        /**
         * Check if this row has identical copy number pattern as another row
         */
        public boolean hasIdenticalPattern(CNVRow other) {
            if (this.copyNumbers.length != other.copyNumbers.length) {
                return false;
            }

            for (int i = 0; i < this.copyNumbers.length; i++) {
                if (this.copyNumbers[i] != other.copyNumbers[i]) {
                    return false;
                }
            }
            return true;
        }

        /**
         * Check if this row can be merged with the next row (same chromosome and consecutive)
         */
        public boolean canMergeWith(CNVRow next) {
            return this.chromosome.equals(next.chromosome) &&
                    this.hasIdenticalPattern(next) &&
                    this.end + 1 == next.start; // Check if positions are consecutive
        }

        @Override
        public String toString() {
            return String.format("%s\t%d\t%d\t%s",
                    chromosome, start, end, Arrays.toString(copyNumbers));
        }
    }

    /**
     * Represents a compressed segment (multiple consecutive rows with same pattern)
     */
    public static class CompressedSegment {
        public String chromosome;
        public long startPosition;  // Start position of first row
        public long endPosition;    // End position of last row
        public int[] copyPattern;   // Copy number pattern

        public CompressedSegment(String chromosome, long startPos, long endPos, int[] copyPattern) {
            this.chromosome = chromosome;
            this.startPosition = startPos;
            this.endPosition = endPos;
            this.copyPattern = copyPattern.clone();
        }

        /**
         * Output in GINKGO format: Chr Start End Cell.1 Cell.2 Cell.3 ...
         */
        public String toOutputString() {
            StringBuilder sb = new StringBuilder();
            sb.append(chromosome).append("\t")
                    .append(startPosition).append("\t")
                    .append(endPosition);

            for (int copyNum : copyPattern) {
                sb.append("\t").append(copyNum);
            }

            return sb.toString();
        }
    }

    /**
     * MAIN COMPRESSION ALGORITHM
     */
    public static List<CompressedSegment> compressConsecutiveRows(List<CNVRow> rows) {
        if (rows.isEmpty()) {
            return new ArrayList<>();
        }

        List<CompressedSegment> compressedSegments = new ArrayList<>();

        int start = 0; // Start from the very first row

        while (start < rows.size()) {
            CNVRow startRow = rows.get(start);
            int current = start;

            // While loop: check if next row has same pattern and same chromosome
            while (current + 1 < rows.size() &&
                    rows.get(current).canMergeWith(rows.get(current + 1))) {
                current++; // Increment - move to next row
            }

            // Now we have start and stop positions
            int stop = current;
            CNVRow endRow = rows.get(stop);

            // Create compressed segment (removed numRowsCompressed parameter)
            CompressedSegment segment = new CompressedSegment(
                    startRow.chromosome,
                    startRow.start,
                    endRow.end,
                    startRow.copyNumbers
            );

            compressedSegments.add(segment);

            // Continue after this segment
            start = stop + 1;
        }

        return compressedSegments;
    }

    /**
     * Read CNV data from file
     */
    public static List<CNVRow> readCNVFile(String filePath) throws IOException {
        List<CNVRow> rows = new ArrayList<>();

        try (BufferedReader reader = new BufferedReader(new FileReader(filePath))) {
            String headerLine = reader.readLine();
            if (headerLine == null) {
                throw new IOException("File is empty");
            }

            String[] headers = headerLine.split("\t");
            if (headers.length < 4) {
                throw new IOException("Invalid format: Need at least Chr, Start, End, Cell columns");
            }

            int numCells = headers.length - 3; // Subtract Chr, Start, End

            String line;
            int lineNumber = 2; // Start from line 2 (after header)

            while ((line = reader.readLine()) != null) {
                line = line.trim();
                if (line.isEmpty()) continue;

                try {
                    String[] parts = line.split("\t");
                    if (parts.length < 3 + numCells) {
                        System.err.println("Warning: Line " + lineNumber + " has insufficient columns, skipping");
                        continue;
                    }

                    String chr = parts[0].trim();
                    long start = Long.parseLong(parts[1].trim());
                    long end = Long.parseLong(parts[2].trim());

                    int[] copyNumbers = new int[numCells];
                    for (int i = 0; i < numCells; i++) {
                        copyNumbers[i] = Integer.parseInt(parts[3 + i].trim());
                    }

                    rows.add(new CNVRow(chr, start, end, copyNumbers));

                } catch (NumberFormatException e) {
                    System.err.println("Warning: Invalid number format at line " + lineNumber + ", skipping");
                }

                lineNumber++;
            }
        }

        return rows;
    }

    /**
     * Write compressed results to file in GINKGO format
     */
    public static void writeCompressedData(List<CompressedSegment> segments,
                                           String outputPath, String[] cellNames) throws IOException {
        try (PrintWriter writer = new PrintWriter(new FileWriter(outputPath))) {
            // Write header in GINKGO format: Chr Start End Cell.1 Cell.2 Cell.3 ...
            writer.print("Chr\tStart\tEnd");
            for (String cellName : cellNames) {
                writer.print("\t" + cellName);
            }
            writer.println();

            // Write compressed segments in GINKGO format
            for (CompressedSegment segment : segments) {
                writer.println(segment.toOutputString());
            }
        }
    }

    /**
     * Extract cell names from header
     */
    public static String[] extractCellNames(String filePath) throws IOException {
        try (BufferedReader reader = new BufferedReader(new FileReader(filePath))) {
            String headerLine = reader.readLine();
            if (headerLine == null) {
                throw new IOException("File is empty");
            }

            String[] headers = headerLine.split("\t");
            String[] cellNames = new String[headers.length - 3];

            for (int i = 3; i < headers.length; i++) {
                cellNames[i - 3] = headers[i];
            }

            return cellNames;
        }
    }

    /**
     * Main method
     */
    public static void main(String[] args) {
        //file path
        String inputFile = "phylonco-lphy/examples/data/copyNumberSim.txt";
        String outputFile = "phylonco-lphy/examples/data/copyNumberSim_compressed.txt";

        // Allow command line arguments to override default paths
        if (args.length >= 1) {
            inputFile = args[0];
        }
        if (args.length >= 2) {
            outputFile = args[1];
        }

        try {
            System.out.println("=== CNV Data Compression Tool ===");
            System.out.println("Input: " + inputFile);

            // Check if input file exists
            File inputFileObj = new File(inputFile);
            if (!inputFileObj.exists()) {
                System.err.println("Error: Input file does not exist: " + inputFile);
                System.exit(1);
            }

            // Read CNV data
            List<CNVRow> originalData = readCNVFile(inputFile);

            if (originalData.isEmpty()) {
                System.err.println("Error: No valid data found in input file");
                System.exit(1);
            }

            // Get cell names and compress data
            String[] cellNames = extractCellNames(inputFile);
            List<CompressedSegment> compressedData = compressConsecutiveRows(originalData);

            // Write results
            writeCompressedData(compressedData, outputFile, cellNames);

            // Show compression statistics
            System.out.println("Compressed " + originalData.size() + " bins → " + compressedData.size() + " segments");

        } catch (IOException e) {
            System.err.println("Error reading/writing files: " + e.getMessage());
            e.printStackTrace();
            System.exit(1);
        } catch (Exception e) {
            System.err.println("Unexpected error: " + e.getMessage());
            e.printStackTrace();
            System.exit(1);
        }
    }
}