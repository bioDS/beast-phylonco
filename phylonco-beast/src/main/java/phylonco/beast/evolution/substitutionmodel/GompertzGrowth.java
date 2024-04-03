package phylonco.beast.evolution.substitutionmodel;

import beast.base.core.BEASTInterface;
import beast.base.core.Description;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.evolution.tree.coalescent.PopulationFunction;
import beast.base.inference.parameter.RealParameter;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.integration.IterativeLegendreGaussIntegrator;

import java.util.ArrayList;
import java.util.List;


    @Description("Coalescent intervals for a gompertz growing population.")
    public class GompertzGrowth extends PopulationFunction.Abstract {
        final public Input<Function> f0Input = new Input<>("f0",
                "Initial proportion of the carrying capacity.", Input.Validate.REQUIRED);
        final public Input<Function> bInput = new Input<>("b",
                "Initial growth rate of tumor growth. Should be greater than 0.", Input.Validate.REQUIRED);
        final public Input<Function> NInfinityInput = new Input<>("NInfinity",
                "Carrying capacity of the population.", Input.Validate.REQUIRED);

        public GompertzGrowth() {
            // Example of setting up inputs with default values
            f0Input.setValue(new RealParameter("0.5"), this);
            bInput.setValue(new RealParameter("0.1"), this);
            NInfinityInput.setValue(new RealParameter("1000"), this);
        }


        private IterativeLegendreGaussIntegrator createIntegrator() {
            int numberOfPoints = 5; // Legendre-Gauss points
            double relativeAccuracy = 1.0e-10; // relative precision
            double absoluteAccuracy = 1.0e-9; // absolute accuracy
            int minimalIterationCount = 2; // Minimum number of iterations
            int maximalIterationCount = 100000; //Maximum number of iterations, adjust as needed
            return new IterativeLegendreGaussIntegrator(numberOfPoints, relativeAccuracy, absoluteAccuracy, minimalIterationCount, maximalIterationCount);
        }



        @Override
        public void initAndValidate() {
            if (f0Input.get() != null && f0Input.get() instanceof RealParameter) {
                RealParameter f0Param = (RealParameter) f0Input.get();
                f0Param.setBounds(Math.max(0.0, f0Param.getLower()), f0Param.getUpper());
            }

            if (bInput.get() != null && bInput.get() instanceof RealParameter) {
                RealParameter bParam = (RealParameter) bInput.get();
                bParam.setBounds(0.0, Double.POSITIVE_INFINITY);  // b should be positive for Gompertz growth
            }

            if (NInfinityInput.get() != null && NInfinityInput.get() instanceof RealParameter) {
                RealParameter NInfinityParam = (RealParameter) NInfinityInput.get();
                NInfinityParam.setBounds(Math.max(0.0, NInfinityParam.getLower()), NInfinityParam.getUpper());
            }

            // Compute N0 from f0 and NInfinity if they are not null and are RealParameters
            if (f0Input.get() != null && NInfinityInput.get() != null &&
                    f0Input.get() instanceof RealParameter && NInfinityInput.get() instanceof RealParameter) {
                double f0 = ((RealParameter) f0Input.get()).getValue();
                double NInfinity = ((RealParameter) NInfinityInput.get()).getValue();
            }
        }

        public double getN0() {
            return getNInfinity() * getF0();
        }

        public double getF0() {
            return f0Input.get().getArrayValue();
        }

        public final double getGrowthRateB() {
            return bInput.get().getArrayValue();
        }

        public double getNInfinity() {
            return NInfinityInput.get().getArrayValue();
        }

        @Override
        public List<String> getParameterIds() {
            List<String> ids = new ArrayList<>();
            // add f0 and b parameter ids
            if (NInfinityInput.get() instanceof BEASTInterface)
                ids.add(((BEASTInterface)NInfinityInput.get()).getID());
            return ids;
        }


        @Override
        public double getPopSize(double t) {
            double f0 = getF0();
            double b = getGrowthRateB();
            double NInfinity = getNInfinity();
            double N0 = getN0();

            return N0 * Math.exp(Math.log(NInfinity / N0) * (1 - Math.exp(b * t)));
        }

        @Override
        public double getIntensity(double t) {
            if (t == 0) return 0;

            UnivariateFunction function = time -> 1 / getPopSize(time);

            //  Use the separate method to create the integrator
            //  return legrandeIntegrator(function, t);
            IterativeLegendreGaussIntegrator integrator = createIntegrator();
            return integrator.integrate(Integer.MAX_VALUE, function, 0, t);
        }

        @Override
        public double getInverseIntensity(double x) {
            return 0;
        }
    }


