package likelihood;

public class TreeLikelihoodWithErrorSlow  extends TreeLikelihoodWithError {
    // slow implementation of sampled error rates

    /**
     * check state for changed variables and update temp results if necessary *
     */
    @Override
    protected boolean requiresRecalculation() {
        if (m_useAmbiguities.get() || m_useTipLikelihoods.get()) {
            setPartials(treeInput.get().getRoot(), dataInput.get().getPatternCount());
        } else {
            setStates(treeInput.get().getRoot(), dataInput.get().getPatternCount());
        }
        return true; // always recalculate
    }

}
