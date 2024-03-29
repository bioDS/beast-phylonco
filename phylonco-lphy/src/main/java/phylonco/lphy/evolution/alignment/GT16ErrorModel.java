package phylonco.lphy.evolution.alignment;

import lphy.base.evolution.alignment.Alignment;
import lphy.base.evolution.alignment.ErrorAlignment;
import lphy.base.evolution.alignment.SimpleAlignment;
import lphy.core.model.GenerativeDistribution;
import lphy.core.model.RandomVariable;
import lphy.core.model.Value;
import lphy.core.model.annotation.Citation;
import lphy.core.model.annotation.GeneratorCategory;
import lphy.core.model.annotation.GeneratorInfo;
import lphy.core.model.annotation.ParameterInfo;
import lphy.core.simulator.RandomUtils;
import org.apache.commons.math3.random.RandomGenerator;
import phylonco.lphy.evolution.datatype.PhasedGenotype;

import java.util.Map;
import java.util.Objects;
import java.util.TreeMap;

/**
 * @author Alexei Drummond
 * @author Kylie Chen
 */
@Citation(
        value = "Alexey Kozlov, Joao Alves, Alexandros Stamatakis, David Posada (2021). CellPhy: accurate and fast probabilistic inference of single-cell phylogenies from scDNA-seq data. bioRxiv 2020.07.31.230292.",
        title = "CellPhy: accurate and fast probabilistic inference of single-cell phylogenies from scDNA-seq data",
        authors = {"Kozlov", "Alves", "Stamatakis", "Posada"},
        year = 2021,
        DOI = "https://doi.org/10.1101/2020.07.31.230292"
)
public class GT16ErrorModel implements GenerativeDistribution<Alignment> {

    Value<Double> epsilon;
    Value<Double> delta;
    Value<Alignment> alignment;

    public final String epsilonParamName = "epsilon";
    public final String deltaParamName = "delta";
    public final String alignmentParamName = "alignment";

    RandomGenerator random;

    // Warning: this line below breaks GUI narrative from displaying correctly!
    // for testing
//    public GT16ErrorModel() {}

    public GT16ErrorModel(
            @ParameterInfo(name = epsilonParamName,
                    narrativeName = "sequencing and amplification error probability",
                    description = "the sequencing and amplification error probability.")
            Value<Double> epsilon,
            @ParameterInfo(name = deltaParamName,
                    narrativeName = "allelic dropout probability",
                    description = "the allelic drop out probability.")
            Value<Double> delta,
            @ParameterInfo(name = alignmentParamName,
                    narrativeName = "genotype alignment",
                    description = "the genotype alignment.")
            Value<Alignment> alignment) {

        if (alignment.value().getSequenceType() != phylonco.lphy.evolution.datatype.PhasedGenotype.INSTANCE) {
            throw new RuntimeException("GT16ErrorModel can only be applied alignments of type " + PhasedGenotype.NAME);
        }

        this.epsilon = epsilon;
        this.delta = delta;
        this.alignment = alignment;
        this.random = RandomUtils.getRandom();
    }

    @Override
    public Map<String, Value> getParams() {
        Map<String, Value> map = new TreeMap<>();
        map.put(epsilonParamName, epsilon);
        map.put(deltaParamName, delta);
        map.put(alignmentParamName, alignment);
        return map;
    }

    @Override
    public void setParam(String paramName, Value value) {
        if (paramName.equals(epsilonParamName)) {
            epsilon = value;
        } else if (paramName.equals(deltaParamName)) {
            delta = value;
        } else if (paramName.equals(alignmentParamName)) {
            alignment = value;
        }
        else throw new RuntimeException("Unrecognised parameter name: " + paramName);
    }

    @GeneratorInfo(
            name = "GT16ErrorModel",
            verbClause = "has",
            narrativeName = "error model",
            category = GeneratorCategory.TAXA_ALIGNMENT,
            description = "An error model distribution on an phased genotype alignment.")
    public RandomVariable<Alignment> sample() {

        Alignment original = alignment.value();
        SimpleAlignment newAlignment = new ErrorAlignment(original.nchar(), original);

        double e = epsilon.value();
        double d = delta.value();

        double[][] errorMatrix = errorMatrix(e, d);

        for (int i = 0; i < newAlignment.ntaxa(); i++) {
            for (int j = 0; j < newAlignment.nchar(); j++) {
                newAlignment.setState(i, j, error(original.getState(i, j), errorMatrix));
            }
        }

        return new RandomVariable<>("D", newAlignment, this);
    }

    public Value<Double> getEpsilon() {
        return getParams().get(epsilonParamName);
    }

    public Value<Double> getDelta() {
        return getParams().get(deltaParamName);
    }

    public Alignment getOriginalAlignment() {
        return Objects.requireNonNull(alignment).value();
    }


    private int error(int state, double[][] errorMatrix) {

        double U = random.nextDouble();

        double[] row = errorMatrix[state];

        double sum = 0;
        for (int i = 0; i < row.length; i++) {
            sum += row[i];
            if (U <= sum) return i;
        }
        throw new RuntimeException("Error in error model! The sum of row should be equal to 1.0");
    }


    public double[][] errorMatrix(double epsilon, double delta) {
        double a = 1 - epsilon + (1/2.0) * delta * epsilon;
        double b = (1 - delta) * (1/6.0) * epsilon;
        double c = (1/6.0) * delta * epsilon;
        double d = (1/2.0) * delta + (1/6.0) * epsilon - (1/3.0) * delta * epsilon;
        double e = (1/6.0) * delta * epsilon;
        double f = (1 - delta) * (1/6.0) * epsilon;
        double g = (1 - delta) * (1 - epsilon);
        // rows are true states Y, columns are observed states X
        // the entries in the matrix are P(X | Y)
        double[][] errorMatrix = {
          // AA AC AG AT CA CC CG CT GA GC GG GT TA TC TG TT
            {a, b, b, b, b, c, 0, 0, b, 0, c, 0, b, 0, 0, c}, // AA
            {d, g, f, f, 0, d, 0, 0, 0, f, e, 0, 0, f, 0, e}, // AC
            {d, f, g, f, 0, e, f, 0, 0, 0, d, 0, 0, 0, f, e}, // AG
            {d, f, f, g, 0, e, 0, f, 0, 0, e, f, 0, 0, 0, d}, // AT
            {d, 0, 0, 0, g, d, f, f, f, 0, e, 0, f, 0, 0, e}, // CA
            {c, b, 0, 0, b, a, b, b, 0, b, c, 0, 0, b, 0, c}, // CC
            {e, 0, f, 0, f, d, g, f, 0, 0, d, 0, 0, 0, f, e}, // CG
            {e, 0, 0, f, f, d, f, g, 0, 0, e, f, 0, 0, 0, d}, // CT
            {d, 0, 0, 0, f, e, 0, 0, g, f, d, f, f, 0, 0, e}, // GA
            {e, f, 0, 0, 0, d, 0, 0, f, g, d, f, 0, f, 0, e}, // GC
            {c, 0, b, 0, 0, c, b, 0, b, b, a, b, 0, 0, b, c}, // GG
            {e, 0, 0, f, 0, e, 0, f, f, f, d, g, 0, 0, 0, d}, // GT
            {d, 0, 0, 0, f, e, 0, 0, f, 0, e, 0, g, f, f, d}, // TA
            {e, f, 0, 0, 0, d, 0, 0, 0, f, e, 0, f, g, f, d}, // TC
            {e, 0, f, 0, 0, e, f, 0, 0, 0, d, 0, f, f, g, d}, // TG
            {c, 0, 0, b, 0, c, 0, b, 0, 0, c, b, b, b, b, a}  // TT
        };
        return errorMatrix;
    }
}