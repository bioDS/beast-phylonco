package phylonco.beast.evolution.populationmodel;


import beast.base.evolution.tree.Tree;
import beast.base.inference.operator.UpDownOperator;

public interface PopFuncWithUpDownOp {
    UpDownOperator getUpDownOperator1(Tree tree);
    UpDownOperator getUpDownOperator2(Tree tree);
}