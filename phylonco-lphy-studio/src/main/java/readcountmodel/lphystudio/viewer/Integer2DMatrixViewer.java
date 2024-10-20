package readcountmodel.lphystudio.viewer;

import lphy.core.model.Value;
import lphystudio.app.graphicalmodelpanel.viewer.Viewer;
import phylonco.lphy.evolution.readcountmodel.Integer2DMatrix;

import javax.swing.*;

public class Integer2DMatrixViewer implements Viewer {
    public Integer2DMatrixViewer() {}
    @Override
    public boolean match(Object value) {
        return value instanceof Integer2DMatrix ||
                (value instanceof Value && ((Value) value).value() instanceof Integer2DMatrix);
    }

    @Override
    public JComponent getViewer(Object value) {
        if (match(value)) {
            return new JTextArea(value.toString());
        }
        String text = ((Value<Integer2DMatrix>) value).value().toString();
        return new JTextArea(text);
    }

    @Override
    public String toString() {
        return "Integer 2D Matrix Viewer";
    }
}
