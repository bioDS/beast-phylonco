//package phylonco.lphybeast.tobeast.values;
//
//import lphy.core.model.Value;
//import lphybeast.BEASTContext;
//import lphybeast.ValueToBEAST;
//import lphybeast.tobeast.values.AlignmentToBEAST;
//import phylonco.lphy.evolution.readcountmodel.MutableAlignment;
//
//
//
//public class MutableAlignmentToBEAST implements ValueToBEAST<MutableAlignment, mutablealignment.MutableAlignment> {
//
//    AlignmentToBEAST alignmentToBEAST = new AlignmentToBEAST();
//
//    @Override
//    public mutablealignment.MutableAlignment valueToBEAST(Value alignmentValue, BEASTContext context) {
//
//        beast.base.evolution.alignment.Alignment beastAlignment = alignmentToBEAST.valueToBEAST(alignmentValue, context);
//
//        if (beastAlignment instanceof mutablealignment.MutableAlignment mutableAlignment) {
//
//            // TODO do your stuff here
//
//            return mutableAlignment;
//        }
//        throw new IllegalArgumentException("The MutableAlignment is required ! But was " + alignmentValue.getType());
//    }
//
//
//    @Override
//    public Class getValueClass() {
//        return MutableAlignment.class;
//    }
//
//    @Override
//    public Class<mutablealignment.MutableAlignment> getBEASTClass() {
//        return mutablealignment.MutableAlignment.class;
//    }
//}
