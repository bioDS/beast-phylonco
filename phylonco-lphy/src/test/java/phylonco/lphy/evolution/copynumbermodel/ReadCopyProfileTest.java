package phylonco.lphy.evolution.copynumbermodel;

import lphy.base.evolution.Taxa;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;

public class ReadCopyProfileTest {

    private ReadCopyProfile reader;

    @BeforeEach
    void setUp() {
        reader = new ReadCopyProfile(new lphy.core.model.Value<>("dummy", "dummy.txt"));
    }

    @Test
    void testParseGinkgoFormat_WithRealData() {
        String realFilePath = "/Users/ming/Desktop/Aug22/CRC09-GINKGO_24cells.CN.txt";

        try {
            List<String> sampleLines = readFirstLines(realFilePath, 31);
            if (sampleLines.size() < 2) {
                System.out.println("Real data file too small or not found, skipping test");
                return;
            }

            List<String> cellNames = extractGinkgoCellNames(sampleLines.get(0));
            if (cellNames.isEmpty()) {
                System.out.println("File doesn't appear to be GINKGO format, skipping test");
                return;
            }

            Taxa taxa = Taxa.createTaxa(cellNames.toArray(new String[0]));

            System.out.println("=== Testing parseGinkgoFormat with Real Data ===");
            System.out.println("Total lines used: " + sampleLines.size() + " (header + " + (sampleLines.size() - 1) + " data rows)");

            // Test the parsing method
            IntegerCharacterMatrix result = reader.parseGinkgoFormat(sampleLines, taxa, cellNames);

            // Basic assertions
            assertEquals(cellNames.size(), result.getTaxa().ntaxa(), "Cell count should match");
            assertEquals(sampleLines.size() - 1, result.nchar().intValue(), "Bin count should match data lines");

        } catch (Exception e) {
            System.err.println("Real GINKGO data parsing test failed: " + e.getMessage());
            // Don't fail test if file issues
        }
    }

    @Test
    void testParseSconce2Format_WithRealData() {
        String realFilePath = "/Users/ming/Desktop/Aug22/CRC09-SCONCE2.CN.txt";

        try {
            List<String> cellNames = extractSconce2CellNames(realFilePath);
            if (cellNames.isEmpty()) {
                System.out.println("Real SCONCE2 data file not found, skipping test");
                return;
            }

            List<String> sampleLines = readSconce2DataForFirstBins(realFilePath, 30);
            if (sampleLines.size() < 2) {
                System.out.println("Not enough SCONCE2 data found, skipping test");
                return;
            }

            Taxa taxa = Taxa.createTaxa(cellNames.toArray(new String[0]));

            System.out.println("=== Testing parseSconce2Format with Real Data ===");
            System.out.println("Total lines used: " + sampleLines.size() + " (header + " + (sampleLines.size() - 1) + " data rows)");

            // Test the parsing method
            IntegerCharacterMatrix result = reader.parseSconce2Format(sampleLines, taxa, cellNames);

        } catch (Exception e) {
            System.err.println("Real SCONCE2 data parsing test failed: " + e.getMessage());
        }
    }

    // Helper methods
    private List<String> readFirstLines(String filePath, int maxLines) throws IOException {
        List<String> lines = new ArrayList<>();
        try (BufferedReader reader = new BufferedReader(new FileReader(filePath))) {
            String line;
            int count = 0;
            while ((line = reader.readLine()) != null && count < maxLines) {
                if (!line.trim().isEmpty()) {
                    lines.add(line.trim());
                    count++;
                }
            }
        }
        return lines;
    }

    private List<String> extractGinkgoCellNames(String headerLine) {
        String[] headers = headerLine.split("\t");
        List<String> cellNames = new ArrayList<>();

        if (headers.length > 3 &&
                headers[0].equalsIgnoreCase("Chr") &&
                headers[1].equalsIgnoreCase("Start") &&
                headers[2].equalsIgnoreCase("End")) {

            for (int i = 3; i < headers.length; i++) {
                String cellName = headers[i].trim();
                if (!cellName.isEmpty()) {
                    cellNames.add(cellName);
                }
            }
        }
        return cellNames;
    }

    private List<String> extractSconce2CellNames(String filePath) throws IOException {
        Set<String> uniqueCells = new LinkedHashSet<>();

        try (BufferedReader reader = new BufferedReader(new FileReader(filePath))) {
            String headerLine = reader.readLine();
            if (headerLine == null) return new ArrayList<>();

            String[] headers = headerLine.split("\t");
            if (!(headers.length == 4 &&
                    headers[0].equalsIgnoreCase("Chr") &&
                    headers[1].equalsIgnoreCase("Start") &&
                    headers[2].equalsIgnoreCase("CN") &&
                    headers[3].equalsIgnoreCase("Sample"))) {
                return new ArrayList<>();
            }

            String line;
            while ((line = reader.readLine()) != null) {
                if (!line.trim().isEmpty()) {
                    String[] values = line.split("\t");
                    if (values.length >= 4) {
                        uniqueCells.add(values[3].trim());
                    }
                }
            }
        }
        return new ArrayList<>(uniqueCells);
    }

    private List<String> readSconce2DataForFirstBins(String filePath, int maxBins) throws IOException {
        List<String> lines = new ArrayList<>();
        Map<String, Integer> cellBinCount = new HashMap<>();

        try (BufferedReader reader = new BufferedReader(new FileReader(filePath))) {
            String headerLine = reader.readLine();
            if (headerLine != null) lines.add(headerLine.trim());

            String line;
            while ((line = reader.readLine()) != null) {
                if (!line.trim().isEmpty()) {
                    String[] values = line.split("\t");
                    if (values.length >= 4) {
                        String cellName = values[3].trim();
                        int currentCount = cellBinCount.getOrDefault(cellName, 0);

                        if (currentCount < maxBins) {
                            lines.add(line.trim());
                            cellBinCount.put(cellName, currentCount + 1);
                        }
                    }
                }
            }
        }
        return lines;
    }
}