//package phylonco.lphybeast.tobeast.generators;
//
//import beast.base.core.BEASTInterface;
//import beast.base.evolution.alignment.Taxon;
//import beast.base.evolution.alignment.TaxonSet;
//import beast.base.inference.parameter.RealParameter;
//import lphy.base.evolution.branchrate.LocalClock;
//import lphy.base.evolution.tree.TimeTree;
//import lphy.base.evolution.tree.TimeTreeNode;
//import lphy.core.model.Value;
//import lphybeast.BEASTContext;
//import lphybeast.GeneratorToBEAST;
//import mf.beast.evolution.branchratemodel.FlexibleLocalClockModel;
//import mf.beast.evolution.branchratemodel.StrictCladeModel;
//import mf.beast.evolution.branchratemodel.StrictLineageClockModel;
//
//import java.util.ArrayList;
//import java.util.List;
//
//public class LocalClockToBeast implements GeneratorToBEAST<LocalClock, FlexibleLocalClockModel> {
//    public FlexibleLocalClockModel generatorToBEAST(LocalClock localClock, BEASTInterface value, BEASTContext context) {
//    /* <branchRateModel spec='mf.beast.evolution.branchratemodel.FlexibleLocalClockModel' id="branchRates"
//                     tree='@Tree.t:fluA'>
//        <rootClockModel id="rate.old"
//                        spec="mf.beast.evolution.branchratemodel.StrictLineageClockModel" clock.rate="@clockRate.old">
//        </rootClockModel>
//
//        <cladeClockModel id="rate.new"
//                         spec="mf.beast.evolution.branchratemodel.StrictCladeModel" taxonset="@taxonset.new"
//                   includeStem="true" clock.rate="@clockRate.new">
//        </cladeClockModel>
//    </branchRateModel>
//    * */
//        FlexibleLocalClockModel flc = new FlexibleLocalClockModel();
//
//        Value<Double> rootRate = localClock.getRootRate();
//        StrictLineageClockModel rootCladeModel = new StrictLineageClockModel();
//        rootCladeModel.initByName("clock.rate", new RealParameter(rootRate.valueToString()));
//
//        Value<Object[]> clades = localClock.getClades();
//        Value<Double[]> cladeRates = localClock.getCladeRates();
//        Value<Boolean> includeStem = localClock.getIncludeStem();
//
//        StrictCladeModel cladeModel = new StrictCladeModel();
//        // TODO: current FLC fix assumes one clade model
//        // add multi clade model in the future
//        Object[] cladesValue = clades.value();
//        for (int i = 0; i < cladesValue.length; i++) {
//            TimeTreeNode node = (TimeTreeNode) cladesValue[i];
//            TaxonSet cladeTaxa = new TaxonSet();
//            ArrayList<Taxon> taxa = new ArrayList<>();
//            List<TimeTreeNode> leaves = node.getAllLeafNodes();
//            for (int j = 0; j < leaves.size(); j++) {
//                Taxon taxon = new Taxon(leaves.get(j).getId());
//                taxa.add(taxon);
//            }
//            cladeTaxa.initByName("taxon", taxa);
//            cladeModel.initByName(
//                    "taxonset", cladeTaxa,
//                    "clock.rate", cladeRates.value()[i].toString(),
//                    "includeStem", includeStem.value()
//            );
//        }
//
//        Value<TimeTree> tree = localClock.getTree();
//
//        flc.initByName(
//                "rootClockModel", rootCladeModel,
//                "cladeClockModel", cladeModel,
//                "tree", context.getBEASTObject(tree).toString()
//        );
//
//
//        flc.initAndValidate();
//
//        return flc;
//    }
//
//    @Override
//    public Class<LocalClock> getGeneratorClass() {
//        return LocalClock.class;
//    }
//
//    @Override
//    public Class<FlexibleLocalClockModel> getBEASTClass() {
//        return FlexibleLocalClockModel.class;
//    }
//}
//
