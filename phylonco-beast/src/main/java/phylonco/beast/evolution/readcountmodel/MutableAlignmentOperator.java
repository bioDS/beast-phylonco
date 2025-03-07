package phylonco.beast.evolution.readcountmodel;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.inference.Operator;
import mutablealignment.MutableAlignment;

import java.util.Random;

@Description("This operator changes an alignment.")
public class MutableAlignmentOperator extends Operator {
    public Input<MutableAlignment> mutableAlignmentInput = new Input<>("mutableAlignment", "mutableAlignment");
    private Random random;
    private MutableAlignment mutableAlignment;
    private Integer site;
    private Integer taxa;
    private Integer oldState;
    private Integer newState;

    public MutableAlignmentOperator() {}


    @Override
    public void initAndValidate() {
        mutableAlignment = mutableAlignmentInput.get();
        random = new Random();
    }


    @Override
    public double proposal() {
        taxa = random.nextInt(mutableAlignment.getTaxonCount());
        site = random.nextInt(mutableAlignment.getSiteCount());
        oldState = mutableAlignment.getSiteValue(taxa, site);
        do {
            newState = random.nextInt(mutableAlignment.getDataType().getStateCount());
        } while (newState == oldState);
        mutableAlignment.setSiteValue(taxa, site, newState);
        return 0;
    }





}
