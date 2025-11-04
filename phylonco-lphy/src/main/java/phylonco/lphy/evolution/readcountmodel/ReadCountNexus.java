package phylonco.lphy.evolution.readcountmodel;

import lphy.core.logger.TextFileFormatted;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;

public class ReadCountNexus implements TextFileFormatted {
    final ReadCountData readCountData;

    public ReadCountNexus(ReadCountData source) {
        readCountData = source;
    }

    @Override
    public int hashCode() {
        return super.hashCode();
    }

    @Override
    public String toString() {
        return super.toString();
    }

    @Override
    public List<String> getTextForFile() {
        // not using this
        return List.of("hello");
    }

    @Override
    public void writeToFile(BufferedWriter writer) {
        try {
            // build taxa
            writer.write(getTaxa(readCountData));
            // build body
            writer.write(getChars(readCountData));
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    public String getTaxa(ReadCountData readCountData) {
        /*
        #NEXUS

        begin taxa;
        dimensions ntax=1;
        taxlabels NC_000001.11;
        end;

         */
        StringBuilder builder = new StringBuilder();
        builder.append("#NEXUS").append("\n\n");

        builder.append("begin taxa;").append("\n");
        builder.append("\tdimensions ntax=").
                append(readCountData.getTaxaNames().length).
                append(";").append("\n");

        builder.append("\ttaxlabels ");
        for (int i = 0; i < readCountData.getTaxaNames().length - 1; i++) {
            builder.append(readCountData.getTaxaNames()[i]).append(" ");
        }
        builder.append(readCountData.getTaxaNames()[readCountData.getTaxaNames().length - 1]).
                append(";").append("\n");
        builder.append("end;").append("\n\n");

        return builder.toString();
    }

    public String getChars(ReadCountData readCountData){
        StringBuilder sb = new StringBuilder();

        sb.append("begin characters;\tdimensions nchar=").
                append(readCountData.getReadCountDataMatrix().length).
                append(";").append("\n");

        sb.append("\tformat datatype=readCount;").append("\n");
        sb.append("\tmatrix").append("\n");

        for (int i = 0; i < readCountData.getTaxaNames().length; i++) {
            sb.append("\t\t").append(readCountData.getTaxaNames()[i]).append(" ").
                    append(writeReadCountLine(readCountData.getReadCountDataMatrix()[i])).append("\n");
        }

        sb.append("end;\n");

        return sb.toString();
    }

    public String writeReadCountLine(ReadCount[] line){
        int n = line.length;
        StringBuilder builder = new StringBuilder();
        for (int i = 0; i < n; i++) {
            int countA = line[i].getCount("A");
            int countC = line[i].getCount("C");
            int countG = line[i].getCount("G");
            int countT = line[i].getCount("T");
            String site = countA + ":" + countC + ":" + countG + ":" + countT;
            builder.append(site);
            if (i < n - 1) {
                builder.append(",");
            }
        }
        builder.append(";");
        return builder.toString();
    }

    @Override
    public String getFileType() {
        return ".nexus";
    }
}
