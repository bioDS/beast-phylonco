package phylonco.beast.evolution.likelihood;

import beast.base.core.Description;
import beast.base.evolution.likelihood.BeerLikelihoodCore;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;

import java.util.List;

@Description("Implementation of optimised tree likelihood calculation with error models")
public class TreeLikelihoodWithErrorFast extends TreeLikelihoodWithError {

    boolean updateLeafPartials = false;

    /**
     * check state for changed variables and update temp results if necessary *
     */
    @Override
    protected boolean requiresRecalculation() {
        hasDirt = Tree.IS_CLEAN;
        if (errorModel != null && errorModel.isDirtyCalculation()) {
            updateLeafPartials = true;
            hasDirt = Tree.IS_FILTHY;
            return true;
        }
        if (dataInput.get().isDirtyCalculation()) {
            hasDirt = Tree.IS_FILTHY;
            return true;
        }
        if (m_siteModel.isDirtyCalculation()) {
            hasDirt = Tree.IS_DIRTY;
            return true;
        }
        if (branchRateModel != null && branchRateModel.isDirtyCalculation()) {
            return true;
        }
        return treeInput.get().somethingIsDirty();
    }

    public void updateLeafPartials() {
        List<Node> leaves = treeInput.get().getExternalNodes();
        for (Node node: leaves) {
            int nodeId = node.getNr();
            likelihoodCore.setNodePartialsForUpdate(nodeId);
            BeerLikelihoodCore beer = (BeerLikelihoodCore) likelihoodCore;
            double[] partials = getLeafPartials(node);
            beer.setCurrentNodePartials(nodeId, partials);
        }
    }

    @Override
    public double calculateLogP() {
        if (updateLeafPartials) {
            updateLeafPartials();
            updateLeafPartials = false;
        }
        return super.calculateLogP();
    }

}
