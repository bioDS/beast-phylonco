package phylonco.beast.evolution.populationmodel;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Loggable;
import beast.base.evolution.tree.coalescent.PopulationFunction;
import beast.base.inference.parameter.RealParameter;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

@Description("Stochastic variable selection for different population growth models.")
public class StochasticVariableSelection extends PopulationFunction.Abstract implements Loggable {
    public final Input<RealParameter> indicatorInput = new Input<>("indicator",
            "The indicator for selecting the population model.", Input.Validate.REQUIRED);
    public final Input<PopulationFunction[]> modelsInput = new Input<>("models",
            "The array of population models.", new PopulationFunction[0]);

    private PopulationFunction selectedModel;

    @Override
    public void initAndValidate() {
        int indicator = (int) indicatorInput.get().getArrayValue();
        PopulationFunction[] models = modelsInput.get();

        if (indicator < 0 || indicator >= models.length) {
            throw new IllegalArgumentException("Invalid indicator value");
        }

        selectedModel = models[indicator];
    }

    @Override
    public double getPopSize(double t) {
        return selectedModel.getPopSize(t);
    }

    @Override
    public double getIntensity(double t) {
        return selectedModel.getIntensity(t);
    }

    //    @Override
    //    public double getIntegral(double start, double finish) {
    //        return selectedModel.getIntegral(start, finish);
    //    }

    @Override
    public double getInverseIntensity(double x) {
        return selectedModel.getInverseIntensity(x);
    }

    @Override
    public List<String> getParameterIds() {
        List<String> ids = new ArrayList<>();
        ids.add(indicatorInput.get().getID());
        for (PopulationFunction model : modelsInput.get()) {
            ids.addAll(model.getParameterIds());
        }
        return ids;
    }

    @Override
    public void init(PrintStream out) {

    }

    @Override
    public void log(long sample, PrintStream out) {

    }

    @Override
    public void close(PrintStream out) {

    }
}
