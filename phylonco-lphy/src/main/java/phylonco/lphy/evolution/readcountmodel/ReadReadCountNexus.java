package phylonco.lphy.evolution.readcountmodel;

import lphy.base.evolution.Taxa;
import lphy.base.evolution.Taxon;
import lphy.base.function.io.ReaderConst;
import lphy.core.io.UserDir;
import lphy.core.logger.LoggerUtils;
import lphy.core.model.DeterministicFunction;
import lphy.core.model.Value;
import lphy.core.model.annotation.GeneratorInfo;
import lphy.core.model.annotation.ParameterInfo;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.NoSuchFileException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;

public class ReadReadCountNexus extends DeterministicFunction<ReadCountData> {

    public ReadReadCountNexus(
            @ParameterInfo(name = ReaderConst.FILE,
                    description = "NEXUS file containing read count data") Value<String> file) {
        if (file == null) {
            throw new IllegalArgumentException("file cannot be null");
        }
        setParam(ReaderConst.FILE, file);
    }

    @GeneratorInfo(name = "readReadCountNexus", examples = {"mpileupToReadCount.lphy"},
            description = "read a nexus file containing read count data")
    @Override
    public Value<ReadCountData> apply() {
        String filePath = getFilePath().value();
        Path path = Paths.get(filePath);

        List<Taxon> taxa = new ArrayList<>();
        List<ReadCount[]> readCountMatrix = new ArrayList<>();
        boolean inMatrix = false;

        try (BufferedReader reader = Files.newBufferedReader(path)) {
            String line;
            while ((line = reader.readLine()) != null) {
                line = line.trim();
                if (line.isEmpty() || line.startsWith("#")) continue;

                if (line.toLowerCase().startsWith("matrix")) {
                    inMatrix = true;
                    continue;
                }
                if (line.toLowerCase().startsWith("end")) {
                    if (inMatrix) break; // stop reading after matrix end
                }

                if (inMatrix) {
                    if (line.endsWith(";"))
                        line = line.substring(0, line.length() - 1);
                    if (line.isEmpty()) continue;

                    String[] parts = line.split("\\s+", 2);
                    if (parts.length < 2) continue;

                    String taxaName = parts[0];
                    Taxon taxon = new Taxon(taxaName);
                    taxa.add(taxon);

                    String[] siteStrings = parts[1].split(",");
                    List<ReadCount> siteList = new ArrayList<>();
                    for (String site : siteStrings) {
                        if (site.isEmpty()) continue;
                        String[] readCount = site.split(":");
                        if (readCount.length != 4) continue;
                        try {
                            int A = Integer.parseInt(readCount[0]);
                            int C = Integer.parseInt(readCount[1]);
                            int G = Integer.parseInt(readCount[2]);
                            int T = Integer.parseInt(readCount[3]);
                            int[] readCounts = new int[]{A,C,G,T};
                            ReadCount rc = new ReadCount(readCounts);
                            siteList.add(rc);
                        } catch (NumberFormatException e) {
                            LoggerUtils.log.warning("Invalid site value for " + taxon + ": " + site);
                        }
                    }
                    readCountMatrix.add(siteList.toArray(new ReadCount[0]));
                }
            }

        } catch (NoSuchFileException e) {
            LoggerUtils.log.severe("File not found: " + e.getMessage() +
                    "\nCurrent working dir = " + UserDir.getUserDir());
        } catch (IOException e) {
            LoggerUtils.logStackTrace(e);
        }

        ReadCountData data = new ReadCountData(
                Taxa.createTaxa(taxa.toArray(new Taxon[0])),
                readCountMatrix.toArray(new ReadCount[0][])
        );

        return new Value<>("", data, this);
    }

    public Value<String> getFilePath() {
        return getParams().get(ReaderConst.FILE);
    }
}
