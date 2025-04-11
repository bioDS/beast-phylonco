package copynumbermodel.lphystudio.viewer;

import lphy.core.model.Value;
import lphystudio.app.graphicalmodelpanel.viewer.Viewer;
import phylonco.lphy.evolution.copynumbermodel.IntegerCharacterMatrix;
import javax.swing.*;

public class TaxaCharacterMatrixViewer implements Viewer {
    public TaxaCharacterMatrixViewer() {}
    @Override
    public boolean match(Object value) {
        return value instanceof IntegerCharacterMatrix ||
                (value instanceof Value && ((Value) value).value() instanceof IntegerCharacterMatrix);
    }

    @Override
    public JComponent getViewer(Object value) {
        if (match(value)) {
            return new JTextArea(value.toString());
        }
        String text = ((Value<IntegerCharacterMatrix>) value).value().toString();
        return new JTextArea(text);
    }

    @Override
    public String toString() {
        return "Integer Matrix Viewer";
    }
}
