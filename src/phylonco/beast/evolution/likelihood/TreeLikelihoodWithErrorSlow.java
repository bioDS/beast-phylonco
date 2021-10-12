package phylonco.beast.evolution.likelihood;

import beast.core.Description;

@Description("Implementation of slow tree likelihood calculation with error models")
public class TreeLikelihoodWithErrorSlow  extends TreeLikelihoodWithError {

    /**
     * always recalculates and updates error matrix
     */
    @Override
    public boolean requiresRecalculation() {
        errorModel.setUpdateFlag(true);
        if (m_useAmbiguities.get() || m_useTipLikelihoods.get()) {
            setPartials(treeInput.get().getRoot(), dataInput.get().getPatternCount());
        } else {
            setStates(treeInput.get().getRoot(), dataInput.get().getPatternCount());
        }
        return true; // always recalculate
    }

}
