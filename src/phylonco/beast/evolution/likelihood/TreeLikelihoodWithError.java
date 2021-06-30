package phylonco.beast.evolution.likelihood;

import beast.core.Description;
import beast.core.Input;
import beast.evolution.alignment.Alignment;
import beast.evolution.likelihood.TreeLikelihood;
import beast.evolution.tree.Node;
import phylonco.beast.evolution.errormodel.ErrorModel;

@Description("Tree likelihood calculation with error models")
public class TreeLikelihoodWithError extends TreeLikelihood {

    final public Input<phylonco.beast.evolution.errormodel.ErrorModel> errorModelInput = new Input<>("errorModel", "error model to use for partials");
    final public Input<Boolean> useTipsEmpiricalInput = new Input<>("useTipsEmpirical", "use tip ambiguities from data", false);

    protected ErrorModel errorModel;
    protected boolean useTipsEmpirical;

    protected boolean useTipLikelihoods = true;
    protected boolean useAmbiguities = true;

    @Override
    public void initAndValidate() {
        // get error model
        errorModel = errorModelInput.get();
        useTipsEmpirical = useTipsEmpiricalInput.get();
        // set fields from TreeLikelihood class
        super.m_useAmbiguities.setValue(useAmbiguities, this);
        super.m_useTipLikelihoods.setValue(useTipLikelihoods, this);
        super.initAndValidate();
    }

    /**
     *
     * @param taxon the taxon name as a string
     * @param data the alignment
     * @return the taxon index of the given taxon name for accessing its sequence data in the given alignment,
     *         or -1 if the taxon is not in the alignment.
     */
    protected int getTaxonIndex(String taxon, Alignment data) {
        int taxonIndex = data.getTaxonIndex(taxon);
        if (taxonIndex == -1) {
            if (taxon.startsWith("'") || taxon.startsWith("\"")) {
                taxonIndex = data.getTaxonIndex(taxon.substring(1, taxon.length() - 1));
            }
            if (taxonIndex == -1) {
                throw new RuntimeException("Could not find sequence " + taxon + " in the alignment");
            }
        }
        return taxonIndex;
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
                    tipLikelihoods = errorModel.getProbabilities(state);
                }
                for (int s = 0; s < nrOfStates; s++) {
                    partials[i] = tipLikelihoods[s];
                    i++;
                }
            }
            likelihoodCore.setNodePartials(node.getNr(), partials);
        } else {
            setPartials(node.getChild(0), nrOfPatterns);
            setPartials(node.getChild(1), nrOfPatterns);
        }
    }

}
