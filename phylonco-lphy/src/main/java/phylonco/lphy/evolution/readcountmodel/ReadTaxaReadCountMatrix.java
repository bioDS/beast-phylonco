package phylonco.lphy.evolution.readcountmodel;


import lphy.base.evolution.Taxa;
import lphy.base.function.io.ReaderConst;
import lphy.core.io.UserDir;
import lphy.core.logger.LoggerUtils;
import lphy.core.model.DeterministicFunction;
import lphy.core.model.Value;
import lphy.core.model.annotation.GeneratorCategory;
import lphy.core.model.annotation.GeneratorInfo;
import lphy.core.model.annotation.IOFunction;
import lphy.core.model.annotation.ParameterInfo;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;

@IOFunction(
        role = IOFunction.Role.dataInput,
        extensions = { ".rc"},
        fileArgument = ReaderConst.FILE
)
public class ReadTaxaReadCountMatrix extends DeterministicFunction<ReadCountData> {

    private Value<Boolean> ref;
    private boolean ifReadReference;
    private int[] refIndex;

    public ReadTaxaReadCountMatrix(@ParameterInfo(name = ReaderConst.FILE, description = "the name of read count matrix file including path, which contains an alignment.") Value<String> filePath,
                                   @ParameterInfo(name = "ref", description = "the optional param(true or false), whether to read reference nucleotide from .rc file. By default, ref is false", optional = true) Value<Boolean> ref
                                   ) {
        if (filePath == null) {
            throw new IllegalArgumentException("The file name can't be null!");
        }
        setParam(ReaderConst.FILE, filePath);
        if (ref != null) {
                this.ref = ref;
        }

    }
    @GeneratorInfo(name="readTaxaReadCountMatrix", verbClause = "is read from", narrativeName = "rc file",
            category = GeneratorCategory.TAXA_ALIGNMENT,
            description = "A function that parses an read count matrix from a rc file.")
    public Value<ReadCountData> apply() {
        ifReadReference = useRef();
        String filePath = ((Value<String>) getParams().get(ReaderConst.FILE)).value();
        Path nexPath = UserDir.getUserPath(filePath);

        try {
            BufferedReader reader = getReader(nexPath.toString());
            ReadCountData data = getReadCountData(reader);
            reader.close();
            return new Value<>(null, data, this);
        } catch (IOException e) {
            throw new RuntimeException("Cannot read file: " + filePath, e);
        }


    }

    private ReadCountData getReadCountData(BufferedReader reader) throws IOException {
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
        Integer[][][] readCounts;
        String[] cellNames;
        int n;
        int l;

        if (ifReadReference) {
            n = allLines.size()-1;
            l = allLines.get(0).split("\t").length-1;
            cellNames = new String[n];
            readCounts = new Integer[n][l][4];
            refIndex = new int[l];
            for (int i = 0; i < l; i++) {
                refIndex[i] = Integer.parseInt(allLines.get(0).split("\t")[i+1]);
            }
            for (int i = 0; i < n; i++) {
                String[] data = allLines.get(i+1).split("\t");
                cellNames[i] = data[0];
                for (int j = 0; j < data.length-1; j++) {
                    String[] readC = data[j+1].split(":");
                    if(readC.length > ReadCount.NUM_NUCLEOTIDES) {
                        throw new IllegalArgumentException("The length of read counts at position: "+j+" is not 4, but is " + readC.length);
                    }
                    for (int k = 0; k < ReadCount.NUM_NUCLEOTIDES; k++) {
                        readCounts[i][j][k] = Integer.parseInt(readC[k]);
                    }
                }
            }
        } else {
            n = allLines.size();
            l = allLines.get(0).split("\t").length-1;
            cellNames = new String[n];
            readCounts = new Integer[n][l][4];
            for (int i = 0; i < n; i++) {
                String[] data = allLines.get(i).split("\t");
                cellNames[i] = data[0];
                for (int j = 0; j < data.length-1; j++) {
                    String[] readC = data[j+1].split(":");
                    if(readC.length != ReadCount.NUM_NUCLEOTIDES) {
                        throw new IllegalArgumentException("The length of read counts at position: "+j+" is not 4, but is " + readC.length);
                    }
                    for (int k = 0; k < ReadCount.NUM_NUCLEOTIDES; k++) {
                        readCounts[i][j][k] = Integer.parseInt(readC[k]);
                    }
                }
            }
        }

        Taxa taxa = Taxa.createTaxa(cellNames);
        ReadCount[][] readCountMatrix = new ReadCount[n][l];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < l; j++) {
                int[] values = new int[ReadCount.NUM_NUCLEOTIDES];
                for (int k = 0; k < ReadCount.NUM_NUCLEOTIDES; k++) {
                    values[k] =  readCounts[i][j][k];
                }
                readCountMatrix[i][j] = new ReadCount(values);

            }
        }
        int[] sitesIndex = new int[l];
        for (int i = 0; i < l; i++) {
            sitesIndex[i] = i;
        }
        ReadCountData readCountData = new ReadCountData(taxa, readCountMatrix, sitesIndex);
        if (ifReadReference) {
            readCountData.setRefIndex(refIndex);
        }
        return readCountData;
    }

    private BufferedReader getReader(String fileName) {
        BufferedReader reader = null;
        try {
            if (!(fileName.endsWith("rc")))
                throw new IOException("Read count matrix file name's suffix is invalid ! " + fileName);

            final Path nexFile = Paths.get(fileName);

            if (!nexFile.toFile().exists() || nexFile.toFile().isDirectory())
                throw new IOException("Cannot find Fasta file ! " + nexFile +
                        ", user.dir = " + System.getProperty("user.dir"));

            reader = Files.newBufferedReader(nexFile); // StandardCharsets.UTF_8
//            reader.mark(READ_AHEAD_LIMIT); // to reset reader back to READ_AHEAD_LIMIT
        } catch (IOException e) {
            LoggerUtils.logStackTrace(e);
        }
        return reader;
    }

    private boolean useRef(){
        if (ref != null) {
            return ref.value();
        } else {
            return false;
        }
    }

}
