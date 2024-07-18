package phylonco.beast.evolution.populationmodel;

import beast.base.core.BEASTInterface;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Loggable;
import beast.base.evolution.tree.coalescent.PopulationFunction;
import beast.base.inference.parameter.IntegerParameter;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

@Description("Stochastic variable selection for different population growth models.")
public class StochasticVariableSelection extends PopulationFunction.Abstract implements Loggable {
    public final Input<IntegerParameter> indicatorInput = new Input<>("indicator",
            "The indicator for selecting the population model.", Input.Validate.REQUIRED);
    public final Input<List<PopulationFunction>> modelsInput = new Input<>("models",
            "The list of population models.", new ArrayList<>());

    private PopulationFunction selectedModel;

    @Override
    public void initAndValidate() {
        if (indicatorInput.get() != null && indicatorInput.get() instanceof IntegerParameter) {
            IntegerParameter IParam = indicatorInput.get();
            IParam.setBounds(Math.max(0, IParam.getLower()), Math.min(3, IParam.getUpper()));
        }

        int indicator = (int) indicatorInput.get().getArrayValue();

        List<PopulationFunction> modelsList = modelsInput.get();

        if (indicator < 0 || indicator >= modelsList.size()) {
            throw new IllegalArgumentException("Invalid indicator value: " + indicator);
        }

        if (modelsList.size() != 4) {
            throw new IllegalArgumentException("There must be exactly 4 population models.");
        }

        selectedModel = modelsList.get(indicator);

        if (selectedModel == null) {
            throw new IllegalArgumentException("Selected model is null. Indicator: " + indicator);
        }
    }

    @Override
    public double getPopSize(double t) {
        return selectedModel.getPopSize(t);
    }

    @Override
    public double getIntensity(double t) {
        return selectedModel.getIntensity(t);
    }

    @Override
    public double getInverseIntensity(double x) {
        return selectedModel.getInverseIntensity(x);
    }

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

    @Override
    public void init(PrintStream out) {}

    @Override
    public void log(long sample, PrintStream out) {}

    @Override
    public void close(PrintStream out) {}
}