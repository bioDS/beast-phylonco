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
    void testParseGinkgoFormat() {
        String filePath = "examples/data/copyNumberSim_ginkgo.txt";

        try {
            List<String> sampleLines = readFirstLines(filePath, 31);
            assertTrue(sampleLines.size() >= 2,
                    "GINKGO files require a minimum of 2 lines (header + data)");

            List<String> cellNames = extractGinkgoCellNames(sampleLines.get(0));

            assertFalse(cellNames.isEmpty(),
                    "GINKGO data file must have valid cell names");

            Taxa taxa = Taxa.createTaxa(cellNames.toArray(new String[0]));

            System.out.println("=== Testing parseGinkgoFormat ===");
            System.out.println("Total lines used: " + sampleLines.size() + " (header + " + (sampleLines.size() - 1) + " data rows)");

            // Test the parsing method
            IntegerCharacterMatrix result = reader.parseGinkgoFormat(sampleLines, taxa, cellNames);

            // Basic assertions
            assertEquals(cellNames.size(), result.getTaxa().ntaxa(), "Cell count should match");
            assertEquals(sampleLines.size() - 1, result.nchar().intValue(), "Bin count should match data lines");

        } catch (IOException e) {
            fail("IOException for file: " + filePath + " - " + e.getMessage());
        }
    }

    @Test
    void testParseSconce2Format() {
        String filePath = "examples/data/copyNumberSim_sconce2.txt";

        try {
            List<String> cellNames = extractSconce2CellNames(filePath);

            assertFalse(cellNames.isEmpty(),
                    "SCONCE2 data file must have valid cell names");

            List<String> sampleLines = readSconce2DataForFirstBins(filePath, 30);

            assertTrue(sampleLines.size() >= 2,
                    "SCONCE2 files require a minimum of 2 lines (header + data)");

            Taxa taxa = Taxa.createTaxa(cellNames.toArray(new String[0]));

            System.out.println("=== Testing parseSconce2Format ===");
            System.out.println("Total lines used: " + sampleLines.size() + " (header + " + (sampleLines.size() - 1) + " data rows)");

            // Test the parsing method
            IntegerCharacterMatrix result = reader.parseSconce2Format(sampleLines, taxa, cellNames);

            assertEquals(cellNames.size(), result.getTaxa().ntaxa(), "Cell count should match");

        } catch (IOException e) {
            fail("IOException for file: " + filePath + " - " + e.getMessage());
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
        // System.out.println("DEBUG: Raw header = '" + headerLine + "'");
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