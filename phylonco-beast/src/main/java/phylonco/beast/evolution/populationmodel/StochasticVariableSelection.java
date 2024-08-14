package phylonco.beast.evolution.populationmodel;

import beast.base.core.BEASTInterface;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Loggable;
import beast.base.evolution.tree.coalescent.PopulationFunction;
import beast.base.inference.StateNode;
import beast.base.inference.parameter.IntegerParameter;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;
import java.util.logging.Logger;

@Description("Stochastic variable selection for different population growth models.")
public class StochasticVariableSelection extends PopulationFunction.Abstract implements Loggable {
    // Logger for logging information and debugging
    private static final Logger logger = Logger.getLogger(StochasticVariableSelection.class.getName());

    // Input for the indicator parameter to select the population model
    public final Input<IntegerParameter> indicatorInput = new Input<>("indicator",
            "The indicator for selecting the population model.", Input.Validate.REQUIRED);
    // Input for the list of population models
    public final Input<List<PopulationFunction>> modelsInput = new Input<>("models",
            "The list of population models.", new ArrayList<>());

    // Selected population model based on the indicator
    private PopulationFunction selectedModel;
    // Keeps track of the previous indicator value to detect changes
    private int previousIndicator = -1;

    // Initialization and validation method
    @Override
    public void initAndValidate() {
        // Ensure the indicator parameter is within bounds
        if (indicatorInput.get() != null && indicatorInput.get() instanceof IntegerParameter) {
            IntegerParameter IParam = indicatorInput.get();
            IParam.setBounds(Math.max(0, IParam.getLower()), Math.min(3, IParam.getUpper()));
        }

        // Get the current indicator value
        selectModel();
    }

    private void selectModel() {
        int indicator = indicatorInput.get().getValue();

        // Get the list of population models
        List<PopulationFunction> modelsList = modelsInput.get();

        // Select the model based on the indicator
        selectedModel = modelsList.get(indicator);

        // Validate the selected model
        if (selectedModel == null) {
            throw new IllegalArgumentException("Selected model is null. Indicator: " + indicator);
        }

        // Log the current and previous indicator values
        //logger.info("Indicator: " + indicator + ", Previous Indicator: " + previousIndicator);

        // Trigger recalculation if the indicator has changed
        if (indicator != previousIndicator) {
            previousIndicator = indicator;
            markLikelihoodRecalculation();
           // logger.info("Model changed, recalculation triggered.");
        }
    }

    // Marks the likelihood for recalculation if the selected model is a StateNode
    private void markLikelihoodRecalculation() {
        if (selectedModel instanceof StateNode) {
            ((StateNode) selectedModel).setEverythingDirty(true);
        }
    }

    // Methods to get population size, intensity, and inverse intensity
    @Override
    public double getPopSize(double t) {
        selectModel();
        return selectedModel.getPopSize(t);
    }

    @Override
    public double getIntensity(double t) {
        selectModel();
        return selectedModel.getIntensity(t);
    }

    @Override
    public double getInverseIntensity(double x) {
        selectModel();
        return selectedModel.getInverseIntensity(x);
    }

    // Method to get the list of parameter IDs
    @Override
    public List<String> getParameterIds() {
        List<String> ids = new ArrayList<>();
        if (indicatorInput.get() instanceof BEASTInterface) {
            ids.add(((BEASTInterface) indicatorInput.get()).getID());
        }
        for (PopulationFunction model : modelsInput.get()) {
            if (model instanceof BEASTInterface) {
                ids.add(((BEASTInterface) model).getID());
            }
        }
        return ids;
    }

    // Loggable interface methods (init, log, close)
    @Override
    public void init(PrintStream out) {}

    @Override
    public void log(long sample, PrintStream out) {}

    @Override
    public void close(PrintStream out) {}

    // Always returns true indicating that the model requires recalculation
    @Override
    public boolean requiresRecalculation() {
        return true;
    }
}
