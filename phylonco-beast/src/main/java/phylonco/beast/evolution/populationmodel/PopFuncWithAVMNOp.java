package phylonco.beast.evolution.populationmodel;

import beast.base.evolution.operator.kernel.AdaptableVarianceMultivariateNormalOperator;
import beast.base.evolution.tree.Tree;

public interface PopFuncWithAVMNOp {

    AdaptableVarianceMultivariateNormalOperator getAVMNOperator(Tree tree);
}
