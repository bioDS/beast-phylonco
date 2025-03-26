//package phylonco.lphy.evolution.readcountmodel;
//
//import lphy.base.evolution.alignment.Alignment;
//import lphy.base.evolution.alignment.SimpleAlignment;
//import lphy.core.model.GenerativeDistribution;
//import lphy.core.model.RandomVariable;
//import lphy.core.model.Value;
//import lphy.core.model.annotation.GeneratorCategory;
//import lphy.core.model.annotation.GeneratorInfo;
//import lphy.core.model.annotation.ParameterInfo;
//
//import java.util.Map;
//import java.util.TreeMap;
//
//public class MutableAlignmentModel implements GenerativeDistribution<Alignment> {
//
//    Value<Alignment> alignment;
//
//    public static final String simpleAlignmetParamName = "A";
//
//
//    public MutableAlignmentModel(@ParameterInfo(name = simpleAlignmetParamName, narrativeName = "simple alignment", description = "the input alingment") Value<Alignment> alignment) {
//
//        this.alignment = alignment;
//    }
//
//    @GeneratorInfo(
//            name = "MutableAlignmentModel",
//            narrativeName = "mutable alignment model",
//            category = GeneratorCategory.TAXA_ALIGNMENT,
//            description = "A model which can transform alignment to mutable alignment."
//    )
//
//    @Override
//    public RandomVariable<Alignment> sample() {
//        Alignment original = alignment.value();
//        SimpleAlignment newAlignment = new MutableAlignment(original.nchar(), original);
//
//
//        for (int i = 0; i < newAlignment.ntaxa(); i++) {
//            for (int j = 0; j < newAlignment.nchar(); j++) {
//                newAlignment.setState(i, j, original.getState(i, j));
//            }
//        }
//
//        return new RandomVariable<>("D", newAlignment, this);
//    }
//
//    @Override
//    public Map<String, Value> getParams() {
//        Map<String, Value> map = new TreeMap<>();
//        map.put(simpleAlignmetParamName, alignment);
//        return map;
//    }
//
//    @Override
//    public void setParam(String paramName, Value value) {
//        if (paramName.equals(simpleAlignmetParamName)) {
//            alignment = value;
//        }
//        else throw new RuntimeException("Unrecognised parameter name: " + paramName);
//    }
//
//    public Value<Alignment> getInputAlignment() {
//        return getParams().get(simpleAlignmetParamName);
//    }
//
//}
