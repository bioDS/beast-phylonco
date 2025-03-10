package phylonco.lphybeast.tobeast.values;

import beast.base.evolution.alignment.Sequence;
import beast.base.evolution.alignment.Taxon;
import beast.base.evolution.alignment.TaxonSet;
import beast.base.evolution.datatype.DataType;
import beast.base.evolution.datatype.UserDataType;
import beast.base.evolution.tree.TraitSet;
import jebl.evolution.sequences.SequenceType;
import lphy.base.evolution.alignment.AlignmentUtils;
import lphy.base.evolution.alignment.SimpleAlignment;
import lphy.base.evolution.datatype.Standard;
import lphy.core.logger.LoggerUtils;
import lphy.core.model.Value;
import lphybeast.BEASTContext;
import lphybeast.ValueToBEAST;
import lphybeast.tobeast.DataTypeUtils;
import phylonco.lphy.evolution.readcountmodel.MutableAlignment;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;


public class MutableAlignmentToBEAST implements ValueToBEAST<MutableAlignment, mutablealignment.MutableAlignment> {

    private static final String DISCRETE = "discrete";

    @Override
    public mutablealignment.MutableAlignment valueToBEAST(Value<MutableAlignment> alignmentValue, BEASTContext context) {

        MutableAlignment alignment = alignmentValue.value();
        SequenceType lphyDataType = alignment.getSequenceType();
        String[] taxaNames = alignment.getTaxaNames();

        mutablealignment.MutableAlignment beastAlignment;
        // TODO BEAST special data types: StandardData, UserDataType, IntegerData
        // 1. Trait Alignment, always 1 site
        if (lphyDataType instanceof Standard standard && alignment.nchar()==1) {
            DataType beastDataType = DataTypeUtils.getUserDataType(standard, true);
            // AlignmentFromTrait
            beastAlignment = new mutablealignment.MutableAlignment();
            // Input<DataType.Base> userDataTypeInput
            beastAlignment.setInputValue("userDataType", beastDataType);

            List<Taxon> taxonList = context.createTaxonList(List.of(taxaNames));
            String traitStr = createTraitString(alignment);

            TraitSet traitSet = new TraitSet();
            traitSet.setInputValue("traitname", DISCRETE);
            traitSet.setInputValue("value", traitStr);

            TaxonSet taxa = new TaxonSet();
            taxa.setInputValue("taxon", taxonList);
            taxa.initAndValidate();

            traitSet.setInputValue("taxa", taxa);
            traitSet.initAndValidate();

            beastAlignment.setInputValue("traitSet", traitSet);
            beastAlignment.initAndValidate();

        } else {
            DataType beastDataType = DataTypeUtils.getBEASTDataType(lphyDataType, context.getDataTypeMap());
            System.out.println("LPhy data type " + lphyDataType + " convert to BEAST data type " + beastDataType);
            // 2. nucleotide, protein, ...
            // sequences
            List<Sequence> sequences = new ArrayList<>();

            int cca = context.getLPhyBeastConfig().compressConstantAlignment;
            // 2.1 trigger get mark[], and check if compress constant sites
            if (cca > 0) {

                Map<Integer, Integer> counter = AlignmentUtils.
                        getConstantSiteWeights(alignment, cca == 2);

                // index is site, value is state, -1 is variable site
//                int[] mark = alignment.getConstantSitesMark();
//                if (alignment.hasConstantSite()) {
//                    long consSiteNum =  Arrays.stream(mark).filter(m -> m != SimpleAlignment.VAR_SITE_STATE).count();
//                    LoggerUtils.log.info("Discover " + consSiteNum + " constant sites.");
//                    if (consSiteNum == alignment.nchar())
//                        throw new RuntimeException("The alignment sites cannot be all constant ! " +
//                                "constant sites " + consSiteNum + " == alignment sites " + alignment.nchar());
//                }
            }

            for (int i = 0; i < taxaNames.length; i++) {
                context.addTaxon(taxaNames[i]);
                // 2.2 if not hasConstantSite(), this returns original sequence
                // otherwise, constant sites are removed
                String s = alignment.getSequenceVarSite(i);
                if (s.length() < 1)
                    throw new RuntimeException("The sequence length cannot < 1 ! Stop at taxon " + taxaNames[i]);
                if ( alignment.hasConstantSite() && i==0 )
                    LoggerUtils.log.info("Keep " + s.length() + " variable sites from the original " + alignment.nchar() + " sites.");
                sequences.add(createBEASTSequence(taxaNames[i], s));
            }
            // normal Alignment
            beastAlignment = new mutablealignment.MutableAlignment();
            beastAlignment.setInputValue("sequence", sequences);

            // 3. morphological data, needs extra <userDataType section
            if (beastDataType instanceof UserDataType) {
                // StandardData.getTypeDescription()
                beastAlignment.setInputValue("dataType", "standard");
                //TODO add FilteredAlignment for standard data ?
                beastAlignment.setInputValue("userDataType", beastDataType);

            } else {
                // Input<String> dataTypeInput
                beastAlignment.setInputValue("dataType", beastDataType.getTypeDescription());
            }

            beastAlignment.initAndValidate();

            // 4 if isCompressConstantSites, then return FilteredAlignment
//            if (context.getLPhyBeastConfig().compressConstantSites && alignment.hasConstantSite()) {
//TODO                // https://www.beast2.org/2019/07/18/ascertainment-correction.html
//                FilteredAlignment filteredAlignment = new FilteredAlignment();
//
//                filteredAlignment.setInputValue("data", beastAlignment);
//                filteredAlignment.setInputValue("filter", "-");
//                // A, C, G, T
//                String weights = createConstantSiteWeights(alignment);
//                filteredAlignment.setInputValue("constantSiteWeights", weights);
//                filteredAlignment.initAndValidate();
//
//                // using LPhy var as ID allows multiple alignments
//                if (!alignmentValue.isAnonymous()) {
//                    beastAlignment.setID("original-" + alignmentValue.getCanonicalId());
//                    filteredAlignment.setID(alignmentValue.getCanonicalId());
//                }
//                // can only use single thread per ascertained alignment
//                return filteredAlignment;
//            }

        }

        // using LPhy var as ID allows multiple alignments
        if (!alignmentValue.isAnonymous()) beastAlignment.setID(alignmentValue.getCanonicalId());

        return beastAlignment;
    }

    // taxa names = traits, ...
    private String createTraitString(SimpleAlignment alignment) {
        StringBuilder builder = new StringBuilder();
        for (int i = 0; i < alignment.ntaxa(); i++) {
            if (i > 0) builder.append(", ");
            builder.append(alignment.getTaxonName(i));
            builder.append("=");
            builder.append(alignment.getSequence(i));
        }
        return builder.toString();
    }

    private Sequence createBEASTSequence(String taxon, String sequence) {
        Sequence seq = new Sequence();
        seq.setInputValue("taxon", taxon);
        seq.setInputValue("value", sequence);
        seq.initAndValidate();
        return seq;
    }

//    private Map<String, Integer> getConstantSiteWeights(SimpleAlignment alignment) {
//
//        String s = alignment.getSequenceVarSite(i);
//
//
//        // index is the state, value is the count
//        NavigableMap<Integer, Integer> counter = AlignmentUtils.countConstantSites(alignment);
//
//        int maxState = counter.lastKey();
//        if (maxState >= alignment.getCanonicalStateCount())
//            LoggerUtils.log.warning("Ambiguous state " + maxState + " is found, count = " + counter.get(maxState));
//
//        List<String> weights = new ArrayList<>();
//        // no ambiguous
//        for (int s = 0; s < alignment.getCanonicalStateCount(); s++) {
//            if (counter.containsKey(s))
//                weights.add(counter.get(s).toString());
//            else
//                weights.add("0");
//        }
//        // check ordering e.g. A, C, G, T
//        SequenceType sequenceType = alignment.getSequenceType();
//        LoggerUtils.log.info("Creating weights mapping to the canonical states in order : " + sequenceType.getCanonicalStates());
//
//        return String.join(" ", weights);
//    }


    @Override
    public Class getValueClass() {
        return MutableAlignment.class;
    }

    @Override
    public Class<mutablealignment.MutableAlignment> getBEASTClass() {
        return mutablealignment.MutableAlignment.class;
    }
}
