package phylonco.lphy.evolution.readcountmodel;

import lphy.base.distribution.DirichletMultinomial;
import lphy.base.evolution.Taxa;
import lphy.base.evolution.alignment.Alignment;
import lphy.core.model.GenerativeDistribution;
import lphy.core.model.RandomVariable;
import lphy.core.model.Value;
import lphy.core.model.annotation.GeneratorInfo;
import lphy.core.model.annotation.ParameterInfo;
import phylonco.lphy.evolution.datatype.PhasedGenotype;
import phylonco.lphy.evolution.datatype.PhasedGenotypeState;

import java.util.Map;
// lphy script
// n = 40;
// l = 200;
// t ~ LogNormal(mean, sd);
// v ~ LogNormal(mean2, sd2);
// ones ~ rep(1, n*l);
// alpha ~ Multinomial(n=1, p=[delta], replicates=n*l) + ones; // alpha_ij = 1 or 2
// sPrime ~ Lognormal(sMean, sStdDev, replicates=n);
// s ~ transpose(rep(sPrime, l));
// mean ~ alpha * rep(t, n*l) * s;
// variance ~
// p ~ // NegativeBinomial p
// r ~ // NegativeBinomial r
// cov ~ NegativeBinomial(...);
// cov ~ reshape(cov, size=[n, l]);
// data ~ PhyloCTMC();
// Reads ~ ReadCountModel(data, cov, delta, epsilon);

public class ReadCountModel implements GenerativeDistribution<ReadCountData> {

    private Value<Integer[][]> coverage;
    private Value<Alignment> data;
    private Value<Double> w;
    private Value<Double> epsilon;
    private RandomVariable<ReadCountData> readCountData;
    private Value<Integer[][]> alpha;
    private DirichletMultinomial dirichletMultinomial = new DirichletMultinomial();

    public static final String dParamName = "D";
    public static final String covParamName = "coverage";
    public static final String alphaParamName = "alpha";
    public static final String epsilonParamName = "epsilon";
    public static final String wParamName = "w";


    public ReadCountModel(
            @ParameterInfo(name = dParamName, description = "genotype alignment.") Value<Alignment> data,
            @ParameterInfo(name = covParamName, description = "coverage.") Value<Integer[][]> coverage,
            @ParameterInfo(name = alphaParamName, description = "allelic dropout events for each cell at each site.") Value<Integer[][]> alpha,
            @ParameterInfo(name = epsilonParamName, description = "sequencing and amplification error probability.") Value<Double> epsilon,
            @ParameterInfo(name = wParamName, description = "overdispersion parameter of Dirichlet multinomial distribution.") Value<Double> w
            ) {
        super();
        this.data = data;
        this.coverage = coverage;
        this.alpha = alpha;
        this.epsilon = epsilon;
        this.w = w;
    }

    @GeneratorInfo(name = "ReadCountModel",
            narrativeName = "read count model",
            description = "A model to simulate the read counts at each site in each cell.")

    // D ~ GT16 genotypes (Alignment)
    // t ~ ...
    // return ACGT {{1,2,3,4}, {4,5,6,7}, ... }

    @Override
    public RandomVariable<ReadCountData> sample() {
        int n = data.value().getTaxa().length();
        int l = data.value().nchar();
        Double eps = epsilon.value();
        Integer[][][] readC = new Integer[n][l][4];
        dirichletMultinomial.setParam("w", this.w);

        Double[][] propensities = {
                {(1 - eps), (eps/3), (eps/3), (eps/3)},             // AA or A_ 0
                {(0.5 - eps/6), (0.5 - eps/6), (eps/6), (eps/6)},   // AC or CA 1
                {(0.5 - eps/6), (eps/6), (0.5 - eps/6), (eps/6)},   // AG or GA 2
                {(0.5 - eps/6),(eps/6),(eps/6),(0.5 - eps/6)},      // AT or TA 3
                {(eps/3), (1 - eps), (eps/3), (eps/3)},             // CC or C_ 4
                {(eps/6), (0.5 - eps/6), (0.5 - eps/6), (eps/6)},   // CG or GC 5
                {(eps/6), (0.5 - eps/6), (eps/6), (0.5 - eps/6)},   // CT or TC 6
                {(eps/3), (eps/3), (1 - eps), (eps/3)},             // GG or G_ 7
                {(eps/6), (eps/6), (0.5 - eps/6), (0.5 - eps/6)},   // GT or TG 8
                {(eps/3), (eps/3), (eps/3), (1 - eps)},             // TT or T_ 9
        };

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < l; j++) {
                int stateIndex = data.value().getState(i, j);
                PhasedGenotypeState genotypeState = PhasedGenotype.getCanonicalState(stateIndex);
                String genotype = genotypeState.getFullName();//get genotype
                Integer cov_ij = coverage.value()[i][j];
                Double[] prob;
                Value<Integer> cover = new Value<>("cover", cov_ij);
                dirichletMultinomial.setParam("n", cover);
                if (genotype.equals("AA")){
                    prob = propensities[0];
                    readC[i][j] = DirichletMultinomialSample(prob);
                }


                else if (genotype.equals("AC") || genotype.equals("CA")){
                    if (alpha.value()[i][j] == 1){
                        if (Math.random() < 0.5){
                            prob = propensities[0];
                            readC[i][j] = DirichletMultinomialSample(prob);
                        }else {
                            prob = propensities[4];
                            readC[i][j] = DirichletMultinomialSample(prob);
                        }
                    }else {
                        prob = propensities[1];
                        readC[i][j] = DirichletMultinomialSample(prob);
                    }


                } else if (genotype.equals("AG") || genotype.equals("GA")){
                    if (alpha.value()[i][j] == 1){
                        if (Math.random() < 0.5){
                            prob = propensities[0];
                            readC[i][j] = DirichletMultinomialSample(prob);
                        }else {
                            prob = propensities[7];
                            readC[i][j] = DirichletMultinomialSample(prob);
                        }
                    }else {
                        prob = propensities[2];
                        readC[i][j] = DirichletMultinomialSample(prob);
                    }


                } else if (genotype.equals("AT") || genotype.equals("TA")){
                    if (alpha.value()[i][j] == 1){
                        if (Math.random() < 0.5){
                            prob = propensities[0];
                            readC[i][j] = DirichletMultinomialSample(prob);
                        }else {
                            prob = propensities[9];
                            readC[i][j] = DirichletMultinomialSample(prob);
                        }
                    }else {
                        prob = propensities[3];
                        readC[i][j] = DirichletMultinomialSample(prob);
                    }


                } else if (genotype.equals("CC")){
                    prob = propensities[4];
                    readC[i][j] = DirichletMultinomialSample(prob);
                }


                else if (genotype.equals("CG") || genotype.equals("GC")){
                    if (alpha.value()[i][j] == 1){
                        if (Math.random() < 0.5){
                            prob = propensities[4];
                            readC[i][j] = DirichletMultinomialSample(prob);
                        }else {
                            prob = propensities[7];
                            readC[i][j] = DirichletMultinomialSample(prob);
                        }
                    } else {
                        prob = propensities[5];
                        readC[i][j] = DirichletMultinomialSample(prob);
                    }


                } else if (genotype.equals("CT") || genotype.equals("TC")){
                    if (alpha.value()[i][j] == 1){
                        if (Math.random() < 0.5){
                            prob = propensities[4];
                            readC[i][j] = DirichletMultinomialSample(prob);
                        }else {
                            prob = propensities[9];
                            readC[i][j] = DirichletMultinomialSample(prob);
                        }
                    } else {
                        prob = propensities[6];
                        readC[i][j] = DirichletMultinomialSample(prob);
                    }


                } else if (genotype.equals("GG")){
                    prob = propensities[7];
                    readC[i][j] = DirichletMultinomialSample(prob);
                }


                else if (genotype.equals("GT") || genotype.equals("TG")){
                    if (alpha.value()[i][j] == 1){
                        if (Math.random() < 0.5){
                            prob = propensities[7];
                            readC[i][j] = DirichletMultinomialSample(prob);
                        }else {
                            prob = propensities[9];
                            readC[i][j] = DirichletMultinomialSample(prob);
                        }
                    }else {
                        prob = propensities[8];
                        readC[i][j] = DirichletMultinomialSample(prob);
                    }


                } else if (genotype.equals("TT")){
                    prob = propensities[9];
                    readC[i][j] = DirichletMultinomialSample(prob);
                }

            }

        }
        Taxa taxa = data.value().getTaxa();
        ReadCount[][] readCountMatrix = new ReadCount[n][l];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < l; j++) {
                int[] values = new int[ReadCount.NUM_NUCLEOTIDES];
                for (int k = 0; k < ReadCount.NUM_NUCLEOTIDES; k++) {
                    values[k] =  readC[i][j][k];
                }
                readCountMatrix[i][j] = new ReadCount(values);

            }
        }
        readCountData = new RandomVariable<>("readCountData", new ReadCountData(taxa, readCountMatrix), this);
        return readCountData;
    }

    private Integer[] DirichletMultinomialSample(Double[] prob){
        Value<Double[]> pro = new Value<>("proMatrix", prob);
        dirichletMultinomial.setParam("p", pro);
        return dirichletMultinomial.sample().value();
    }

    @Override
    public Map<String, Value> getParams() {
        return Map.of(
                dParamName, data,
                covParamName, coverage,
                alphaParamName, alpha,
                epsilonParamName, epsilon,
                wParamName, w
        );
    }

    @Override
    public void setParam(String paramName, Value value) {
        if (paramName.equals(dParamName)) data = value;
        else if (paramName.equals(covParamName)) coverage = value;
        else if (paramName.equals(alphaParamName)) alpha = value;
        else if (paramName.equals(epsilonParamName)) epsilon = value;
        else if (paramName.equals(wParamName)) w = value;
        else throw new RuntimeException("Unrecognised parameter name: " + paramName);

        //super.setParam(paramName, value); // constructDistribution
    }
}
