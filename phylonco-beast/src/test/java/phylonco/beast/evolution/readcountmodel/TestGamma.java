package phylonco.beast.evolution.readcountmodel;

import org.apache.commons.math3.special.Gamma;

public class TestGamma {

    public static void testOne() {
        double x = 142;
        double gamma = Gamma.gamma(x);
        double logGamma = Gamma.logGamma(x);
        double gamma2 = Math.exp(logGamma);

        System.out.println("gamma(" + x + ") = " + gamma);
        System.out.println("logGamma(" + x + ") = " + logGamma);
        System.out.println("gamma2 = exp( logGamma(" + x + ") ) = " + gamma2);
        double y = Gamma.gamma(141) * 141;
        System.out.println("y = gamma(141) * 141 = " + y);

        System.out.println("log( gamma(" + y + ") ) = " + Math.log(y));
//        System.out.println("exp( logGamma(" + y + ") ) = " + Math.exp(Gamma.logGamma(y)));

        // https://www.wolframalpha.com/input?i=gamma%28142%29
        // Input: Γ(142)
        // Result (decimal approximation): 1.8981437590761709694285264141107677937281750118953493737972... × 10^243

    }

    public static void testTwo() {

        double x = 300;
        double gamma = Gamma.gamma(x);
        double logGamma = Gamma.logGamma(x);
        double gamma2 = Math.exp(logGamma);

        System.out.println("gamma(" + x + ") = " + gamma);
        System.out.println("logGamma(" + x + ") = " + logGamma);
        System.out.println("gamma2 = exp( logGamma(" + x + ") ) = " + gamma2);
//        double y = Gamma.gamma(141) * 141;
//        System.out.println("y = gamma(141) * 141 = " + y);
//
//        System.out.println("log( gamma(" + y + ") ) = " + Math.log(y));
    }



    public static void main(String[] args) {

        testTwo();

    }
}
