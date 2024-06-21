package phylonco.beast.evolution.populationmodel;



public class GompertzToLogistic {
    public static double computeT50(double NInfinity, double N0, double b) {
        if (N0 >= NInfinity || b <= 0) {
            throw new IllegalArgumentException("N0 must be less than NInfinity and b must be greater than 0.");
        }
        double ratio = NInfinity / N0;
        double proportion = 0.5;
        double t50 = Math.log(1 - Math.log(proportion) / Math.log(ratio)) / b;
        return t50;
    }

    public static void main(String[] args) {
        double f0Min = 0.9033885;
        double f0Max = 0.9977269;
        double bMin = 0.0962412;
        double bMax = 0.1361323;
        double step = 0.01;
        double N0 = 2000000000000000.0;

        double minT50 = Double.MAX_VALUE;
        double maxT50 = Double.MIN_VALUE;

        for (double f0 = f0Min; f0 <= f0Max; f0 += step) {
            double NInfinity = N0 / f0;
            for (double b = bMin; b <= bMax; b += step) {
                double t50 = computeT50(NInfinity, N0, b);
                if (t50 < minT50) {
                    minT50 = t50;
                }
                if (t50 > maxT50) {
                    maxT50 = t50;
                }
                System.out.printf("f0: %.2f, b: %.2f -> t50: %.3f\n", f0, b, t50);
            }
        }

        // 检查边界值
        double[] f0Edges = {f0Min, f0Max};
        double[] bEdges = {bMin, bMax};
        for (double f0 : f0Edges) {
            double NInfinity = N0 / f0;
            for (double b : bEdges) {
                double t50 = computeT50(NInfinity, N0, b);
                if (t50 < minT50) {
                    minT50 = t50;
                }
                if (t50 > maxT50) {
                    maxT50 = t50;
                }
                System.out.printf("Boundary - f0: %.2f, b: %.2f -> t50: %.3f\n", f0, b, t50);
            }
        }

        System.out.printf("t50 range: %.3f to %.3f\n", minT50, maxT50);
    }
}
