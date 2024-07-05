package phylonco.beast.evolution.populationmodel;

import beast.base.core.Description;
import beast.base.inference.Operator;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.util.Randomizer;
import lphy.core.model.annotation.GeneratorInfo;
import lphy.core.model.annotation.ParameterInfo;

/**
 * Stochastic variable selection operator for MCMC sampling using Indicator
 */
@Description("Stochastic variable selection operator for MCMC sampling using Indicator")
public class StochasticVariableSelectionOperator extends Operator {

    private final ConstantGrowth constantGrowth;
    private final ExponentialGrowth exponentialGrowth;
    private final LogisticGrowth logisticGrowth;
    private final GompertzGrowth_t50 gompertzGrowth;
    private final IntegerParameter indicator;

    public StochasticVariableSelectionOperator(
            @ParameterInfo(name = "constantGrowth", description = "Constant growth model") ConstantGrowth constantGrowth,
            @ParameterInfo(name = "exponentialGrowth", description = "Exponential growth model") ExponentialGrowth exponentialGrowth,
            @ParameterInfo(name = "logisticGrowth", description = "Logistic growth model") LogisticGrowth logisticGrowth,
            @ParameterInfo(name = "gompertzGrowth", description = "Gompertz growth model") GompertzGrowth_t50 gompertzGrowth,
            @ParameterInfo(name = "Indicator", description = "Indicator for model selection") IntegerParameter indicator) {
        this.constantGrowth = constantGrowth;
        this.exponentialGrowth = exponentialGrowth;
        this.logisticGrowth = logisticGrowth;
        this.gompertzGrowth = gompertzGrowth;
        this.indicator = indicator;
    }

    private static final int WINDOW_SIZE = 1;

    @Override
    public void initAndValidate() {
        // No initialization needed for this operator
    }

    @Override
    public double proposal() {
        int currentValue = indicator.getValue();
        int newValue = proposeNewIndicator(indicator);

        // Update the indicator value
        indicator.setValue(newValue);

        return 0.0;
    }

    /**
     * Proposes a new value for the indicator within the defined window size.
     *
     * @param indicator The IntegerParameter representing the indicator.
     * @return The new proposed value for the indicator.
     */
    private int proposeNewIndicator(IntegerParameter indicator) {
        int currentValue = indicator.getValue();
        int newValue = currentValue + Randomizer.nextInt(2 * WINDOW_SIZE + 1) - WINDOW_SIZE;

        // Ensure the new value is within valid bounds
        if (newValue < indicator.getLower() || newValue > indicator.getUpper()) {
            return currentValue; // Retain the original value if out of bounds
        }
        return newValue;
    }

    /**
     * SVSPopfunction gets the population size at time t based on the current indicator value.
     *
     * @param modelIndex The index of the selected model.
     * @param model The selected model.
     * @param t The time at which to get the population size.
     * @return The population size at time t.
     */
    @GeneratorInfo(name="SVSPopfunction",
            narrativeName = "Kingman's coalescent tree prior",
            description="The Kingman coalescent distribution over tip-labelled time trees.")
    public double SVSPopfunction(int modelIndex, Object model, double t) {
        switch (modelIndex) {
            case 0:
                return ((ConstantGrowth) model).getPopSize(t);
            case 1:
                return ((ExponentialGrowth) model).getPopSize(t);
            case 2:
                return ((LogisticGrowth) model).getPopSize(t);
            case 3:
                return ((GompertzGrowth_t50) model).getPopSize(t);
            default:
                throw new IllegalArgumentException("Invalid modelIndex value");
        }
    }


    }

