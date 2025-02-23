package phylonco.beast.evolution.populationmodel;

import beast.base.evolution.tree.Tree;
import beast.base.inference.operator.UpDownOperator;

public interface PopFuncWithUpOp {
    UpDownOperator getUpOperator(Tree tree);
}