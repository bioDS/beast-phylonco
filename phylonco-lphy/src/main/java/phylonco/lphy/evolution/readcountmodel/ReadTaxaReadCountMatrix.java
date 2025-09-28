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
import java.io.FileReader;
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

    public ReadTaxaReadCountMatrix(@ParameterInfo(name = ReaderConst.FILE, description = "the name of read count matrix file including path, which contains an alignment.") Value<String> filePath) {
        if (filePath == null) throw new IllegalArgumentException("The file name can't be null!");
        setParam(ReaderConst.FILE, filePath);

    }
    @GeneratorInfo(name="readTaxaReadCountMatrix", verbClause = "is read from", narrativeName = "rc file",
            category = GeneratorCategory.TAXA_ALIGNMENT,
            description = "A function that parses an read count matrix from a rc file.")
    public Value<ReadCountData> apply() {

        String filePath = ((Value<String>) getParams().get(ReaderConst.FILE)).value();
        Path nexPath = UserDir.getUserPath(filePath);

        try {
            BufferedReader reader = new BufferedReader(new FileReader(nexPath.toFile()));
            ReadCountData data = getReadCountData(reader);
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
        int n = allLines.size();
        int l = allLines.get(0).split("\t").length-1;
        String[] cellNames = new String[n];
        Integer[][][] readCounts = new Integer[n][l][4];
        for (int i = 0; i < n; i++) {
            String[] data = allLines.get(i).split("\t");
            cellNames[i] = data[0];
            for (int j = 0; j < data.length-1; j++) {
                String[] readC = data[j+1].split(":");
                if(readC.length > ReadCount.NUM_NUCLEOTIDES) {System.out.println("the number of read counts at position: "+j+" for cell: "+ data[0] +" is higher than 4");}
                for (int k = 0; k < ReadCount.NUM_NUCLEOTIDES; k++) {
                    readCounts[i][j][k] = Integer.parseInt(readC[k]);
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
        return readCountData;
    }

    private BufferedReader getReader(String fileName) {
        BufferedReader reader = null;
        try {
            // Common Extensions: .fasta, .fas, .fa
            // Other Extensions (depending on sequence type):
            // .fna (nucleotide sequences), .ffn (nucleotide sequences),
            // .faa (amino acid sequences), .frn (RNA sequences), .mpfa (multiple protein FASTA)
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

}
