package phylonco.lphy.evolution.alignment;

import jebl.evolution.sequences.SequenceType;
import lphy.base.evolution.Taxa;
import lphy.base.evolution.alignment.Alignment;
import lphy.base.evolution.alignment.SimpleAlignment;
import lphy.core.logger.LoggerUtils;
import lphy.core.logger.TextFileFormatted;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

public class HaploidAlignment extends SimpleAlignment implements TextFileFormatted {

    public HaploidAlignment(Map<String, Integer> idMap, int nchar, SequenceType sequenceType) {
        super(idMap, nchar, sequenceType);
    }

    public HaploidAlignment(Taxa taxa, int nchar, SequenceType sequenceType) {
        super(taxa, nchar, sequenceType);
    }

    public HaploidAlignment(int nchar, Alignment source) {
        super(nchar, source);
    }

    @Override
    public List<String> getTextForFile() {
        List<String> lines = new ArrayList<>();
        List<String> Haploid1 = new ArrayList<>();
        List<String> Haploid2 = new ArrayList<>();

        for (int i=0; i < ntaxa(); i++) {
            try {
                String taxonName = getTaxaNames()[i];
                Haploid1.add(">" + taxonName + " | the first haploid sequence");

                String sequence = getSequence(i);
                StringBuilder even = new StringBuilder();
                StringBuilder odd = new StringBuilder();

                for (int s=0;s<sequence.length();s++){
                    if (s % 2 == 0){
                        even.append(sequence.charAt(s));
                    } else {
                        odd.append(sequence.charAt(s));
                    }
                }
                lines.add(even.toString());

                Haploid2.add(">" + taxonName + " | the second haploid sequence");
                lines.add(odd.toString());
            } catch (Exception ex) {
                LoggerUtils.log.severe("Error at " + i + " taxa (" + getTaxaNames()[i] + ") in " +
                        this.getClass().getName());
            }

        }
        return lines;
    }

    @Override
    public String getFileType() {
        return ".fa";
    }
}
