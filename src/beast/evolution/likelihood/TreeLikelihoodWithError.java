package beast.evolution.likelihood;

import beast.core.Input;
import beast.evolution.alignment.Alignment;
import beast.evolution.datatype.DataTypeWithError;
import beast.evolution.likelihood.TreeLikelihood;
import beast.evolution.tree.Node;

public class TreeLikelihoodWithError extends TreeLikelihood {

    final public Input<DataTypeWithError> errorModelInput = new Input<>("errorModel", "error model to use for partials");
    final public Input<Boolean> useTipsEmpiricalInput = new Input<>("useTipsEmpirical", "use tip ambiguities from data", false);

    protected DataTypeWithError errorModel;
    protected boolean useTipsEmpirical; // set to true to use tips from data, otherwise use tips from model

    @Override
    public void initAndValidate() {
        // get error model
        errorModel = errorModelInput.get();
        useTipsEmpirical = useTipsEmpiricalInput.get();
        // set fields from TreeLikelihood class
        super.m_useAmbiguities.set(true);
        super.m_useTipLikelihoods.set(true);
        super.initAndValidate();
    }

    @Override
    protected void setPartials(Node node, int nrOfPatterns) {
        if (node.isLeaf()) {
            Alignment data = dataInput.get();
            int nrOfStates = data.getDataType().getStateCount();
            double[] partials = new double[nrOfPatterns * nrOfStates];
            int t = getTaxonIndex(node.getID(), data); // taxon index
            int i = 0;
            for (int p = 0; p < nrOfPatterns; p++) {
                int state = data.getPattern(t, p);
                double[] tipLikelihoods;
                if (useTipsEmpirical) {
                    tipLikelihoods = data.getTipLikelihoods(t, p);
                } else {
                    tipLikelihoods = errorModelInput.get().getProbabilities(state);
                }
                for (int s = 0; s < nrOfStates; s++) {
                    if (tipLikelihoods == null) {
                        boolean[] stateSet = data.getStateSet(state);
                        if (stateSet[s]) {
                            partials[i] = 1.0;
                        } else {
                            partials[i] = 0.0;
                        }
                    } else {
                        partials[i] = tipLikelihoods[s];
                    }
                    i++;
                }
            }
            likelihoodCore.setNodePartials(node.getNr(), partials);
//            boolean debug = true;
//            if (debug) {
//                System.out.println("partials " + data.getTaxaNames().get(t));
//                int count = 1;
//                for (double p : partials) {
//                    System.out.print(p + " ");
//                    if (count % nrOfStates == 0) {
//                        System.out.println("; ");
//                    }
//                    count++;
//                }
//                System.out.println();
//            }
        } else {
            setPartials(node.getLeft(), nrOfPatterns);
            setPartials(node.getRight(), nrOfPatterns);
        }
    }

}
