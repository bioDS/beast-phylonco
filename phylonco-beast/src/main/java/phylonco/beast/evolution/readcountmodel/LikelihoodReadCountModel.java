package phylonco.beast.evolution.readcountmodel;


import beast.base.core.Input;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.tree.Node;
import beast.base.inference.Distribution;
import beast.base.inference.State;
import beast.base.inference.parameter.RealParameter;
import org.apache.commons.math3.special.Gamma;
import phylonco.beast.evolution.datatype.ReadCount;

import java.util.List;
import java.util.Random;


public class LikelihoodReadCountModel extends Distribution {

    public Input<Alignment> alignmentInput = new Input<>("alignment", "alignment");
    public Input<ReadCount> readCountInput = new Input<>("readCount", "nucleotide read counts");

    // epsilon, allelic dropout, ... parameters
    public Input<RealParameter> epsilonInput = new Input<>("epsilon", "sequencing error");
    public Input<RealParameter> deltaInput = new Input<>("delta", "allelic dropout probability");
    public Input<RealParameter> tInput = new Input<>("t", "mean of allelic coverage");
    public Input<RealParameter> vInput = new Input<>("v", "variance of allelic coverage");
    public Input<RealParameter> sInput = new Input<>("s", "size factor of cell");
    public Input<RealParameter> wInput = new Input<>("w", "overdispersion parameter of Dirichlet multinomial distribution");

    // other parameters of read count model

    private RealParameter epsilon;
    private RealParameter delta;
    private RealParameter t;
    private RealParameter v;
    private RealParameter s;
    private RealParameter w;

    private double alpha1;
    private double alpha2;
    private double mean1;
    private double mean2;
    private double variance1;
    private double variance2;
    private double negp1;
    private double negp2;
    private double negr1;
    private double negr2;
    private Double[][] propensities;



    @Override
    public List<String> getArguments() {
        return List.of();
    }

    @Override
    public List<String> getConditions() {
        return List.of();
    }

    @Override
    public void sample(State state, Random random) {

    }

    @Override
    public void initAndValidate() {
        // checking parameters correct

        // get parameters
        this.epsilon = epsilonInput.get();
        this.delta = deltaInput.get();
        this.t = tInput.get();
        this.v = vInput.get();
        this.s = sInput.get();
        this.w = wInput.get();
        this.alpha1 = 1;
        this.alpha2 = 2;




        calculateLogP();
    }

    public void calculateLogPLeaf(Node node, int[] states) {

    }

    //Calculate the log likelihood of read count model by summarizing the log likelihood at each site
    @Override
    public double calculateLogP() {
        Double eps = epsilon.getValue();
        Double tv = this.t.getValue();
        Double vv = this.v.getValue();
        Double wv = this.w.getValue();
        this.propensities = new Double[][]{
                {(1 - eps), (eps/3), (eps/3), (eps/3)},             // AA or A_ 0
                {(0.5 - eps/6), (0.5 - eps/6), (eps/6), (eps/6)},   // AC or CA 1
                {(0.5 - eps/6), (eps/6), (0.5 - eps/6), (eps/6)},   // AG or GA 2
                {(0.5 - eps/6), (eps/6), (eps/6), (0.5 - eps/6)},   // AT or TA 3
                {(eps/3), (1 - eps), (eps/3), (eps/3)},             // CC or C_ 4
                {(eps/6), (0.5 - eps/6), (0.5 - eps/6), (eps/6)},   // CG or GC 5
                {(eps/6), (0.5 - eps/6), (eps/6), (0.5 - eps/6)},   // CT or TC 6
                {(eps/3), (eps/3), (1 - eps), (eps/3)},             // GG or G_ 7
                {(eps/6), (eps/6), (0.5 - eps/6), (0.5 - eps/6)},   // GT or TG 8
                {(eps/3), (eps/3), (eps/3), (1 - eps)},             // TT or T_ 9
        };
        this.logP = 0;
        for (int i = 0; i < alignmentInput.get().getTaxonCount(); i++) {
            Double sv = this.s.getValues()[i];
            this.mean1 = alpha1 * tv * sv;
            this.mean2 = alpha2 * tv * sv;
            this.variance1 = mean1 + Math.pow(this.alpha1, 2) * vv * Math.pow(sv, 2);
            this.variance2 = mean2 + Math.pow(this.alpha2, 2) * vv * Math.pow(sv, 2);
            this.negp1 = this.mean1 / this.variance1;
            this.negp2 = this.mean2 / this.variance2;
            this.negr1 = Math.pow(this.mean1, 2) / (this.variance1 - this.mean1);
            this.negr2 = Math.pow(this.mean2, 2) / (this.variance2 - this.mean2);
            for (int j = 0; j < alignmentInput.get().getSiteCount(); j++) {///？
                // dirichlet multinomial pmf
                int patternIndex = alignmentInput.get().getPatternIndex(j);
                int genotypeState = alignmentInput.get().getPattern(i, patternIndex);
                int[] readCountNumbers = readCountInput.get().getReadCounts(i, j);
                this.logP += logLiklihoodRC(genotypeState, readCountNumbers, wv);
            }
        }
        //System.out.println("logLikelihood = " + logP + "; t = " + this.t + "; v = " + this.v + "; s = " + this.s + "; w = " + this.w);
        //System.out.println("this.logP: " + this.logP);
        return this.logP;
    }


    // calculate probability of read counts given genotype
    // genotypeState represents genotype alignment
    public double logLiklihoodRC(int genotypeState, int[] readCountNumbers, double w) {
        int coverage = 0;
        for (int i = 0; i < readCountNumbers.length; i++) {
            coverage = coverage + readCountNumbers[i];
        }
        double deltav = delta.getValue();

        int[] indices = getGenotypeIndices(genotypeState);

        double logLikelihood;
        double logLikelihoodDirichletMDDiploid;
        double logLikelihoodDirichletMDHaploid0;
        double logLikelihoodDirichletMDHaploid1;
        double logCoverageLikelihoodDiploid;
        double logCoverageLikelihoodHaploid;
        double part0;
        double part1;
        double part2;
        double max;


        if (homozygous(genotypeState)) {
            logLikelihoodDirichletMDDiploid = logLikelihoodDirichletMD(w, coverage, propensities[indices[0]], readCountNumbers);
            logCoverageLikelihoodDiploid = logCoverageLikelihood(readCountNumbers, negp2, negr2);
            logLikelihoodDirichletMDHaploid0 = logLikelihoodDirichletMD(w, coverage, propensities[indices[1]], readCountNumbers);
            logCoverageLikelihoodHaploid = logCoverageLikelihood(readCountNumbers, negp1, negr1);
            part0 = logLikelihoodDirichletMDDiploid + logCoverageLikelihoodDiploid + Math.log(1 - deltav);
            part1 = logLikelihoodDirichletMDHaploid0 + logCoverageLikelihoodHaploid + Math.log(deltav);
            max = Math.max(part0, part1);
            if (part0 == max){
                logLikelihood = part0 + Math.log(1 + Math.exp(part1 - part0));
            }else {
                logLikelihood = part1 + Math.log(1 + Math.exp(part0 - part1));
            }
        } else {
            logLikelihoodDirichletMDDiploid = logLikelihoodDirichletMD(w, coverage, propensities[indices[0]], readCountNumbers);
            logCoverageLikelihoodDiploid = logCoverageLikelihood(readCountNumbers, negp2, negr2);
            logLikelihoodDirichletMDHaploid0 = logLikelihoodDirichletMD(w, coverage, propensities[indices[1]], readCountNumbers);
            logCoverageLikelihoodHaploid = logCoverageLikelihood(readCountNumbers, negp1, negr1);
            logLikelihoodDirichletMDHaploid1 = logLikelihoodDirichletMD(w, coverage, propensities[indices[2]], readCountNumbers);
            part0 = logLikelihoodDirichletMDDiploid + logCoverageLikelihoodDiploid + Math.log(1 - deltav);
            part1 = Math.log(0.5) + logLikelihoodDirichletMDHaploid0 + logCoverageLikelihoodHaploid + Math.log(deltav);
            part2 = Math.log(0.5) + logLikelihoodDirichletMDHaploid1 + logCoverageLikelihoodHaploid + Math.log(deltav);
            max = Math.max(part0, Math.max(part1, part2));
            if (part0 == max){
                logLikelihood = part0 + Math.log(1 + Math.exp(part1 - part0)) + Math.log(1 + Math.exp(part2 - part0));
            }else if (part1 == max){
                logLikelihood = part1 + Math.log(1 + Math.exp(part0 - part1)) + Math.log(1 + Math.exp(part2 - part1));
            }else {
                logLikelihood = part2 + Math.log(1 + Math.exp(part0 - part2)) + Math.log(1 + Math.exp(part1 - part2));
            }
        }
        return logLikelihood;
    }
    //calculate the probability at each site given read count(coverage)(negative-binomial distribution)
    public double logCoverageLikelihood(int[] readCountNumbers, double p, double r) {
        // negative binomial pmf
        double logCoverageLikelihood;
        int c = 0;
        for (int i = 0; i < readCountNumbers.length; i++) {
            c = c + readCountNumbers[i];
        }
        logCoverageLikelihood = Gamma.logGamma(c + r) - Gamma.logGamma(r) - Gamma.logGamma(c + 1) + r * Math.log(p) + c * Math.log(1 - p);
        return logCoverageLikelihood;
    }
    //get indices from propensities matrix by given genotype
    private static int[] getGenotypeIndices(int genotypeState) {
        int[][] indexTable = {
                {0,0},      //AA (AA and A_)
                {1,0,4},    //AC (AC, A_ and C_)
                {2,0,7},    //AG (AG, A_ and G_)
                {3,0,9},    //AT (AT, A_ and T_)
                {1,4,0},    //CA (CA, A_ and C_)
                {4,4},      //CC (CC and C_)
                {5,4,7},    //CG (CG, C_ and G_)
                {6,4,9},    //CT (CT, C_ and T_)
                {2,7,0},    //GA (GA, G_ and A_)
                {5,7,4},    //GC (GC, G_ and C_)
                {7,7},      //GG (GG and G_)
                {8,7,9},    //GT (GT, G_ and T_)
                {3,9,0},    //TA (TA, T_ and A_)
                {6,9,4},    //TC (TC, T_ and C_)
                {8,9,7},    //TG (TG, T_ and G_)
                {9,9},      //TT (TT and T_)
        };

        int[] indices = indexTable[genotypeState];
        return indices;
    }

    //Determining whether a genotype is homozygous or not
    private boolean homozygous(int genotype){
        return switch (genotype){
            case 0, 5, 10, 15 -> true;
            default -> false;
        };
    }

    //calculate the likelihood given read count (multinomial distribution)
    public double logLikelihoodDirichletMD(double w, int coverage, Double[] propensities, int[] readCountNumbers){
        double logLikelihood = logFFunction(coverage, w);
        for (int i = 0; i < readCountNumbers.length; i++) {
            logLikelihood = logLikelihood - logFFunction(readCountNumbers[i], w * propensities[i]);
        }
        return logLikelihood;
    }

    public double logFFunction(int coverage, double w){
        double result;
        if (coverage > 0){
            result = Math.log(coverage) + Gamma.logGamma(w) + Gamma.logGamma(coverage) - Gamma.logGamma(w + coverage);
            return result;
        } else return 0.0;
    }


}
