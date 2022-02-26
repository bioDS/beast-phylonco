package phylonco.beast.evolution.operators;

import beast.core.Description;
import beast.core.Input;
import beast.core.Operator;
import beast.core.StateNode;
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;

import java.util.*;

@Description("Swap many pairs of parameter values within a single proposal")
public class BlockPairSwapOperator extends Operator {

    final public Input<RealParameter> parameterInput = new Input<>("parameter", "a real parameter to swap individual values for");
    final public Input<IntegerParameter> pairInput = new Input<>("pairs", "a vector specifying which pairs to swap", Input.Validate.REQUIRED);

    RealParameter parameter;
    IntegerParameter pairs;

    int numPairs = 0;

    static Map<Integer, List<Integer>> staticMap;

    @Override
    public void initAndValidate() {
        parameter = parameterInput.get();
        pairs = pairInput.get();

        assert(parameter != null);
        assert(pairs != null);

        pairs.initAndValidate();

        assert(pairs.getDimension() == parameter.getDimension());

        HashMap<Integer, List<Integer>> indexMap = new HashMap<>();
        for (int i = 0; i < parameter.getDimension(); i++) {
            int key = pairs.getValue(i);
            if (key > 0) {
                if (indexMap.containsKey(key)) {
                    List<Integer> indexList = indexMap.get(key);
                    indexList.add(i);
                } else {
                    List<Integer> indexList = new ArrayList();
                    indexList.add(i);
                    indexMap.put(key, indexList);
                }
            }
        }

        // checking pairs have correct sizes
        for (Integer key: indexMap.keySet()) {
            assert(indexMap.get(key).size() == 2);
        }

        staticMap = Collections.unmodifiableMap(indexMap);
        numPairs = staticMap.size();
    }

    @Override
    public double proposal() {
        for (Integer key: staticMap.keySet()) {
            int left = staticMap.get(key).get(0);
            int right = staticMap.get(key).get(1);
            parameter.swap(left, right);
        }

        return 0; // symmetric move return log Hastings ratio
    }

    @Override
    public List<StateNode> listStateNodes() {
        final List<StateNode> list = new ArrayList<>();
        list.add(parameter);
        return list;
    }

}
