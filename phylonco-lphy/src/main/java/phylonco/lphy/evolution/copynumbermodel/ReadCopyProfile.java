package phylonco.lphy.evolution.copynumbermodel;

import lphy.base.evolution.Taxa;
import lphy.base.function.io.ReaderConst;
import lphy.core.io.UserDir;
import lphy.core.model.DeterministicFunction;
import lphy.core.model.Value;
import lphy.core.model.annotation.GeneratorCategory;
import lphy.core.model.annotation.GeneratorInfo;
import lphy.core.model.annotation.IOFunction;
import lphy.core.model.annotation.ParameterInfo;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Path;
import java.util.List;
import java.util.*;

/**
 * D = readCopyProfile(file="realCNVData.txt");
 */
@IOFunction(
        role = IOFunction.Role.dataInput,
        extensions = {".txt"},
        fileArgument = ReaderConst.FILE
)
public class ReadCopyProfile extends DeterministicFunction<IntegerCharacterMatrix> {

    public ReadCopyProfile(@ParameterInfo(name = ReaderConst.FILE, description = "the name of txt file of copy number profiles.") Value<String> filePath) {

        if (filePath == null) {
            throw new IllegalArgumentException("The file path can't be null!");
        }

        setParam(ReaderConst.FILE, filePath);

    }

    @GeneratorInfo(name = "readCopyProfile", verbClause = "is read from", narrativeName = "Copy number profile file",
            category = GeneratorCategory.TAXA_ALIGNMENT,
            description = "A function that parses copy number profiles from txt files (GINKGO, SCONCE2 formats).")
    public Value<IntegerCharacterMatrix> apply() {
        String filePath = ((Value<String>) getParams().get(ReaderConst.FILE)).value();
        Path readPath = UserDir.getUserPath(filePath);

        try {
            BufferedReader reader = new BufferedReader(new FileReader(readPath.toFile()));
            IntegerCharacterMatrix data = getData(reader);
            return new Value<>(null, data, this);
        } catch (IOException e) {
            throw new RuntimeException("Cannot read file: " + filePath, e);
        }

    }

    private IntegerCharacterMatrix getData(BufferedReader reader) throws IOException {
        // TODO: implement this method to parse data from BufferedReader input
        // Step 1: Read all lines from the file
        List<String> allLines = new ArrayList<>();
        String line;

        while ((line = reader.readLine()) != null) {
            if (!line.trim().isEmpty()) { // Skip empty lines
                allLines.add(line.trim());
            }
        }
        if (allLines.isEmpty()) {
            throw new IllegalArgumentException("File is empty or contains no valid data");
        }
        // Step 2: Parse header to determine format
        String headerLine = allLines.get(0);
        String[] headers = headerLine.split("\t");

        // Step 3: Detect format and extract cell names
        List<String> cellNames = new ArrayList<>();
        boolean isGinkgoFormat = false;

        // Check if GINKGO format: Chr, Start, End, Cell.1, Cell.2, ...
        if (headers.length > 3 &&
                headers[0].equalsIgnoreCase("Chr") &&
                headers[1].equalsIgnoreCase("Start") &&
                headers[2].equalsIgnoreCase("End")) {

            isGinkgoFormat = true;
            System.out.println("Detected GINKGO format");

            // Extract cell names from columns 3 onwards
            for (int i = 3; i < headers.length; i++) {
                cellNames.add(headers[i].trim());
            }
        }
        // Check if SCONCE2 format: Chr, Start, CN, Sample
        else if (headers.length == 4 &&
                headers[0].equalsIgnoreCase("Chr") &&
                headers[1].equalsIgnoreCase("Start") &&
                headers[2].equalsIgnoreCase("CN") &&
                headers[3].equalsIgnoreCase("Sample")) {

            System.out.println("Detected SCONCE2 format");
            // Extract unique cell names from Sample column
            Set<String> uniqueCells = new LinkedHashSet<>();
            for (int i = 1; i < allLines.size(); i++) { // Skip header
                String[] values = allLines.get(i).split("\t");
                if (values.length >= 4) {
                    uniqueCells.add(values[3].trim());
                }
            }
            cellNames.addAll(uniqueCells);
        } else {
            throw new IllegalArgumentException("Unknown file format. Expected GINKGO or SCONCE2 format. Found headers: " + Arrays.toString(headers));
        }
        if (cellNames.isEmpty()) {
            throw new IllegalArgumentException("No cell names found in file");
        }
        // System.out.println("Found " + cellNames.size() + " cells: " + cellNames); //Test

        // Create Taxa object
        String[] cellArray = cellNames.toArray(new String[0]);
        Taxa taxa = Taxa.createTaxa(cellArray);

        if (isGinkgoFormat) {
            return parseGinkgoFormat(allLines, taxa, cellNames);
        } else {
            return parseSconce2Format(allLines, taxa, cellNames);
        }
    }


    IntegerCharacterMatrix parseGinkgoFormat(List<String> lines, Taxa taxa, List<String> cellNames) {
        // Skip header line
        List<String> dataLines = lines.subList(1, lines.size());
        int nbins = dataLines.size();
        int numberOfCells = cellNames.size();

        System.out.println("Parsing GINKGO: " + nbins + " bins, " + numberOfCells + " cells"); //Test

        // Create matrix: rows = cells, columns = bins
        IntegerCharacterMatrix matrix = new IntegerCharacterMatrix(taxa, nbins);
        // Process each data line (each line = one genomic bin)
        for (int binIndex = 0; binIndex < nbins; binIndex++) {
            String[] values = dataLines.get(binIndex).split("\t");

            // Check we have enough columns: Chr + Start + End + copy numbers for each cell
            if (values.length < 3 + numberOfCells) {
                throw new IllegalArgumentException("Line " + (binIndex + 2) + " has insufficient columns. Expected " +
                        (3 + numberOfCells) + ", found " + values.length);
            }

            // Extract copy numbers for each cell (skip Chr, Start, End columns)
            for (int cellIndex = 0; cellIndex < numberOfCells; cellIndex++) {
                String copyNumberStr = values[3 + cellIndex].trim();
                try {
                    int copyNumber = Integer.parseInt(copyNumberStr);
                    // Set state: cellIndex (which cell), binIndex (which bin), copyNumber (value)
                    matrix.setState(cellIndex, binIndex, copyNumber);
                } catch (NumberFormatException e) {
                    throw new IllegalArgumentException("Invalid copy number '" + copyNumberStr +
                            "' at line " + (binIndex + 2) +
                            ", cell " + cellNames.get(cellIndex));
                }
            }
        }
        System.out.println("Successfully parsed GINKGO format"); //Test
        return matrix;

    }

    IntegerCharacterMatrix parseSconce2Format(List<String> lines, Taxa taxa, List<String> cellNames) {
        // Step 1: Group data by genomic position (Chr:Start)
        Map<String, List<String>> dataByBin = new LinkedHashMap<>();

        // Skip header and group by genomic position
        for (int i = 1; i < lines.size(); i++) {
            String[] values = lines.get(i).split("\t");
            if (values.length >= 4) {
                String binKey = values[0] + ":" + values[1]; // "Chr:Start"
                dataByBin.computeIfAbsent(binKey, k -> new ArrayList<>()).add(lines.get(i));
            }
        }

        int nbins = dataByBin.size();
        int numberOfCells = cellNames.size();

        System.out.println("Parsing SCONCE2: " + nbins + " bins, " + numberOfCells + " cells");  //Test

        // Create matrix
        IntegerCharacterMatrix matrix = new IntegerCharacterMatrix(taxa, nbins);

        // Step 2: Process each genomic bin
        int binIndex = 0;
        for (Map.Entry<String, List<String>> binEntry : dataByBin.entrySet()) {
            String genomicPosition = binEntry.getKey();
            List<String> binData = binEntry.getValue();

            // Create map: cellName -> copyNumber for this bin
            Map<String, Integer> cellCopyNumbers = new HashMap<>();

            // Parse all rows for this genomic position
            for (String dataLine : binData) {
                String[] values = dataLine.split("\t");
                if (values.length >= 4) {
                    String cellName = values[3].trim();   // Sample column
                    String copyNumberStr = values[2].trim(); // CN column

                    try {
                        int copyNumber = Integer.parseInt(copyNumberStr);
                        cellCopyNumbers.put(cellName, copyNumber);
                    } catch (NumberFormatException e) {
                        throw new IllegalArgumentException("Invalid copy number '" + copyNumberStr +
                                "' for cell " + cellName +
                                " at position " + genomicPosition);
                    }
                }
            }

            // Step 3: Set copy numbers for each cell in this bin
            for (int taxonIndex = 0; taxonIndex < numberOfCells; taxonIndex++) {
                String cellName = cellNames.get(taxonIndex);
                Integer copyNumber = cellCopyNumbers.get(cellName);

                if (copyNumber == null) {
                    throw new IllegalArgumentException("Missing copy number data for cell '" + cellName +
                            "' at genomic position " + genomicPosition);
                }

                // Set state: taxonIndex (which cell), binIndex (which bin), copyNumber (value)
                matrix.setState(taxonIndex, binIndex, copyNumber);
            }

            binIndex++;
        }
        System.out.println("Successfully parsed SCONCE2 format"); //Test
        return matrix;
    }
}

