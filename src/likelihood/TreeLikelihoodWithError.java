package likelihood;

import beast.core.Input;
import beast.evolution.alignment.Alignment;
import beast.evolution.datatype.DataTypeWithError;
import beast.evolution.likelihood.TreeLikelihood;
import beast.evolution.tree.Node;

public class TreeLikelihoodWithError extends TreeLikelihood {
    // set fields from TreeLikelihood class
    // "useAmbiguities": true
    // "useTipLikelihoods": true

    final public Input<DataTypeWithError> errorModelInput = new Input<>("errorModel", "error model to use for partials");
    final public Input<Boolean> useTipsEmpiricalInput = new Input<>("useTipsEmpirical", "use tip ambiguities from data", false);

    protected DataTypeWithError errorModel;
    protected boolean useTipsEmpirical; // set to true to use tips from data, otherwise use tips from model

    @Override
    public void initAndValidate() {
        super.initAndValidate();
        errorModel = errorModelInput.get();
        useTipsEmpirical = useTipsEmpiricalInput.get();
    }

    @Override
    protected void setPartials(Node node, int patternCount) {
        if (node.isLeaf()) {
            Alignment data = dataInput.get();
            int states = data.getDataType().getStateCount();
            double[] partials = new double[patternCount * states];
            int k = 0;
            int taxonIndex = getTaxonIndex(node.getID(), data);
            for (int patternIndex_ = 0; patternIndex_ < patternCount; patternIndex_++) {
                double[] tipLikelihoods;
                if (useTipsEmpirical) {
                    tipLikelihoods = data.getTipLikelihoods(taxonIndex, patternIndex_);
                } else {
                    int stateIndex = data.getPattern(taxonIndex, patternIndex_);
                    tipLikelihoods = errorModelInput.get().getProbabilities(stateIndex);
                }
                if (tipLikelihoods != null) {
                    for (int state = 0; state < states; state++) {
                        partials[k++] = tipLikelihoods[state];
                    }
                } else {
                    int stateIndex = data.getPattern(taxonIndex, patternIndex_);
                    boolean[] stateSet = data.getStateSet(stateIndex);
                    for (int state = 0; state < states; state++) {
                        partials[k++] = (stateSet[state] ? 1.0 : 0.0);
                    }
                }
            }
            likelihoodCore.setNodePartials(node.getNr(), partials);
            // print partials for debugging
            boolean debug = true;
            if (debug) {
                System.out.println("partials " + data.getTaxaNames().get(taxonIndex));
                int count = 1;
                for (double p : partials) {
                    System.out.print(p + " ");
                    if (count % states == 0) {
                        System.out.println("; ");
                    }
                    count++;
                }
                System.out.println();
            }
            // end debug
        } else {
            setPartials(node.getLeft(), patternCount);
            setPartials(node.getRight(), patternCount);
        }
    }

}
