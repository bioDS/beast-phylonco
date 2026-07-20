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
import phylonco.lphy.evolution.datatype.UnphasedGenotype;

import java.util.Map;
import java.util.Objects;
import java.util.TreeMap;

/**
 * @author Kylie Chen
 * GT10 error model from Kozlov et al., 2021.
 */
@Citation(
        value = "Alexey Kozlov, Joao Alves, Alexandros Stamatakis, David Posada (2021). CellPhy: accurate and fast probabilistic inference of single-cell phylogenies from scDNA-seq data. bioRxiv 2020.07.31.230292.",
        title = "CellPhy: accurate and fast probabilistic inference of single-cell phylogenies from scDNA-seq data",
        authors = {"Kozlov", "Alves", "Stamatakis", "Posada"},
        year = 2021,
        DOI = "https://doi.org/10.1101/2020.07.31.230292"
)
public class GT10ErrorModel implements GenerativeDistribution<Alignment> {

    Value<Double> epsilon;
    Value<Double> delta;
    Value<Alignment> alignment;

    public final String epsilonParamName = "epsilon";
    public final String deltaParamName = "delta";
    public final String alignmentParamName = "alignment";

    RandomGenerator random;

    public GT10ErrorModel(
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

        if (alignment.value().getSequenceType() != phylonco.lphy.evolution.datatype.UnphasedGenotype.INSTANCE) {
            throw new RuntimeException("GT10ErrorModel can only be applied alignments of type " + UnphasedGenotype.NAME);
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
            name = "GT10ErrorModel",
            verbClause = "has",
            narrativeName = "error model",
            category = GeneratorCategory.TAXA_ALIGNMENT,
            description = "An error model distribution on an unphased genotype alignment.")
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
        double a = 1 - epsilon + (1/2.0) * delta * epsilon; // scenario I
        double b = (1 - delta) * (1/3.0) * epsilon; // scenario II updated
        double c = (1/6.0) * delta * epsilon; // scenario III
        double d = (1/2.0) * delta + (1/6.0) * epsilon - (1/3.0) * delta * epsilon; // scenario IV
        double e = (1/6.0) * delta * epsilon; // scenario V
        double f = (1 - delta) * (1/6.0) * epsilon; // scenario VI updated
        double g = (1 - delta) * (1 - epsilon); // scenario VII
        // rows are true states Y, columns are observed states X
        // the entries in the matrix are P(X | Y)
        double X = -1;
        double[][] errorMatrix = {
                // AA AC AG AT CC CG CT GG GT TT observed states
                {a, b, b, b, c, 0, 0, c, 0, c}, // AA 0
                {d, g, f, f, d, f, f, e, 0, e}, // AC 1
                {d, f, g, f, e, f, 0, d, f, e}, // AG 2
                // AA AC AG AT CC CG CT GG GT TT
                {d, f, f, g, e, 0, f, e, f, d}, // AT 3
                {c, b, 0, 0, a, b, b, c, 0, c}, // CC 4
                {e, f, f, 0, d, g, f, d, f, e}, // CG 5
                // AA AC AG AT CC CG CT GG GT TT
                {e, f, 0, f, d, f, g, e, f, d}, // CT 6
                {c, 0, b, 0, c, b, 0, a, b, c}, // GG 7
                {e, 0, f, f, e, f, f, d, g, d}, // GT 8
                {c, 0, 0, b, c, 0, b, c, b, a}  // TT 9
        };
        return errorMatrix;
    }
}