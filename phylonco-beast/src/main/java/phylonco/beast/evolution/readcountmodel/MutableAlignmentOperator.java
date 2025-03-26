//package phylonco.beast.evolution.readcountmodel;
//
//import beast.base.core.Description;
//import beast.base.core.Input;
//import beast.base.inference.Operator;
//import beast.base.util.Randomizer;
//import mutablealignment.MutableAlignment;
//
//@Description("This operator changes an alignment.")
//public class MutableAlignmentOperator extends Operator {
//
//
//    public Input<MutableAlignment> mutableAlignmentInput = new Input<>("mutableAlignment", "mutableAlignment");
//    private MutableAlignment mutableAlignment;
//    private int numStates;
//    private int numSites;
//    private int numTaxa;
//
//    public MutableAlignmentOperator() {}
//
//    @Override
//    public void initAndValidate() {
//        mutableAlignment = mutableAlignmentInput.get();
//        numStates = mutableAlignment.getDataType().getStateCount();
//        numSites = mutableAlignment.getSiteCount();
//        numTaxa = mutableAlignment.getTaxonCount();
//    }
//
//
//    @Override
//    public double proposal() {
//        int taxa = Randomizer.nextInt(numTaxa);
//        int site = Randomizer.nextInt(numSites);
//        int oldState = mutableAlignment.getSiteValue(taxa, site);
//        int newState;
//        do {
//            newState = Randomizer.nextInt(numStates);
//        } while (newState==oldState);
//        mutableAlignment.setSiteValue(taxa, site, newState);
//        return 0;
//    }
//
//
//
//}
