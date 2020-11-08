package beast.evolution.likelihood;

// slow implementation of sampled error rates
public class TreeLikelihoodWithErrorSlow  extends TreeLikelihoodWithError {

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
