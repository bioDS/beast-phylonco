package phylonco.beast.evolution.readcountmodel;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.inference.Operator;
import beast.base.util.Randomizer;
import mutablealignment.MutableAlignment;

@Description("This operator changes an alignment.")
public class MutableAlignmentOperator extends Operator {
    public Input<MutableAlignment> mutableAlignmentInput = new Input<>("mutableAlignment", "mutableAlignment");
    private MutableAlignment mutableAlignment;
    private Integer site;
    private Integer taxa;
    private Integer oldState;
    private Integer newState;

    public MutableAlignmentOperator() {}


    @Override
    public void initAndValidate() {
        mutableAlignment = mutableAlignmentInput.get();
    }


    @Override
    public double proposal() {
        taxa = Randomizer.nextInt(mutableAlignment.getTaxonCount());
        site = Randomizer.nextInt(mutableAlignment.getSiteCount());
        oldState = mutableAlignment.getSiteValue(taxa, site);
        do {
            int numStates = mutableAlignment.getDataType().getStateCount();
            newState = Randomizer.nextInt(numStates);
        } while (newState == oldState);
        mutableAlignment.setSiteValue(taxa, site, newState);
        return 0;
    }





}
