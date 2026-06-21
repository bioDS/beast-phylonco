package phylonco.lphybeast.loggerhelper;

import beast.base.core.BEASTInterface;
import beast.base.core.Loggable;
import beast.base.evolution.branchratemodel.BranchRateModel;
import beast.base.evolution.tree.Tree;
import beast.base.inference.Logger;
import beast.base.spec.evolution.sitemodel.SiteModel;
import com.google.common.collect.Multimap;
import lphy.core.model.GraphicalModelNode;
import lphybeast.BEASTContext;
import lphybeast.tobeast.loggers.LoggerHelper;
import phylonco.beast.evolution.datatype.ReadCount;
import phylonco.beast.evolution.likelihood.ReadCountTreeLikelihood;
import phylonco.beast.evolution.logger.SampledGenotypeLogger;
import phylonco.beast.evolution.readcountmodel.LikelihoodReadCountModel;

import java.util.List;

import static java.lang.Math.toIntExact;

public class AlignmentLoggerHelper implements LoggerHelper {
    final protected BEASTContext context;
    final protected ReadCountTreeLikelihood treeLikelihood;
    final protected SiteModel siteModel;
    final protected LikelihoodReadCountModel readCountModel;
    final protected ReadCount readCount;
    private SampledGenotypeLogger sampledGenotypeLogger;
    String fileName;

    public AlignmentLoggerHelper(ReadCountTreeLikelihood treeLikelihood, SiteModel siteModel, LikelihoodReadCountModel readCountModel,
                                 ReadCount readCount, BEASTContext context) {
        this.context = context;
        this.treeLikelihood = treeLikelihood;
        this.siteModel = siteModel;
        this.readCountModel = readCountModel;
        this.readCount = readCount;
    }


    @Override
    public Logger createLogger(long l, Multimap<BEASTInterface, GraphicalModelNode<?>> multimap) {
        sampledGenotypeLogger = new SampledGenotypeLogger();
        sampledGenotypeLogger.setInputValue("tree", (Tree) treeLikelihood.treeInput.get());
        sampledGenotypeLogger.setInputValue("siteModel", siteModel);
        BranchRateModel branchRateModel = treeLikelihood.branchRateModelInput.get();
        if (branchRateModel != null) {
            sampledGenotypeLogger.setInputValue("branchRateModel", branchRateModel);
        }
        sampledGenotypeLogger.setInputValue("readCountModel", readCountModel);
        sampledGenotypeLogger.setInputValue("readCount", readCount);
        sampledGenotypeLogger.initAndValidate();
        sampledGenotypeLogger.setID("sampledGenotype");
        Logger logger = new Logger();
        // Must convert to int
        logger.setInputValue("logEvery", toIntExact(l));
        logger.setInputValue("log", sampledGenotypeLogger);
        logger.setInputValue("fileName", getFileName());
        return logger;
    }

    @Override
    public List<Loggable> getLoggables() {
        List<Loggable> loggables = List.of();
        loggables.add(sampledGenotypeLogger);
        return loggables;
    }

    @Override
    public String getFileName() {
        return fileName;
    }

    @Override
    public void setFileName(String fileStem, boolean isMultiple) {
            fileName = fileStem + "_alignment.log";
    }
}
