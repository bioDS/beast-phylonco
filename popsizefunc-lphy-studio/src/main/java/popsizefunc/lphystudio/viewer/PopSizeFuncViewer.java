package popsizefunc.lphystudio.viewer;

import lphy.core.model.Value;
import lphystudio.app.graphicalmodelpanel.viewer.Viewer;
import popsizefunc.lphy.evolution.popsize.PopulationFunction;

import javax.swing.*;

public class PopSizeFuncViewer implements Viewer {

    /**
     * Required by ServiceLoader.
     */
    public PopSizeFuncViewer() {
    }

    @Override
    public boolean match(Object value) {
        return value instanceof PopulationFunction ||
                (value instanceof Value && ((Value) value).value() instanceof PopulationFunction);
    }

    @Override
    public JComponent getViewer(Object value) {
        //TODO create PopulationFunctionComponent((Value<PopulationFunction>) value)
        if (value instanceof PopulationFunction) {
            return new JLabel("PopulationFunction TODO : " + value.toString());
        }
        return new JLabel("Value<PopulationFunction> TODO : " + ((Value<PopulationFunction>) value).value().toString());
    }
    @Override
    public String toString() {
        return "Pop-size function Viewer";
    }

}
