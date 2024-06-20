package readcountmodel.lphystudio.viewer;

import lphy.base.evolution.coalescent.PopulationFunction;
import lphy.core.model.Value;
import lphystudio.app.graphicalmodelpanel.viewer.Viewer;
import phylonco.lphy.evolution.readcountmodel.ReadCountData;
import phylonco.lphy.evolution.readcountmodel.ReadCountModel;

import javax.swing.*;

public class ReadCountModelViewer implements Viewer{
    public ReadCountModelViewer(){}

    @Override
    public boolean match(Object value) {
        return value instanceof ReadCountData ||
                (value instanceof Value && ((Value) value).value() instanceof ReadCountData);

    }

    @Override
    public JComponent getViewer(Object value) {
        if (match(value)) {
            return new JTextArea(value.toString());
        }
        String text = ((Value<ReadCountData>) value).value().toString();
        return new JTextArea(text);
    }

    @Override
    public String toString() { return "Read Count Viewer"; }
}

