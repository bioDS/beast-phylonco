package phylonco.lphybeast.tobeast.values;

import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.alignment.Sequence;
import lphy.core.model.Value;
import lphybeast.BEASTContext;
import lphybeast.ValueToBEAST;
import phylonco.lphy.evolution.copynumbermodel.*;

import java.util.ArrayList;
import java.util.List;

public class IntegerCharacterMatrixToBEAST implements ValueToBEAST<IntegerCharacterMatrix, Alignment> {

    @Override
    public Alignment valueToBEAST(Value<IntegerCharacterMatrix> value, BEASTContext context) {
        // Extract the IntegerCharacterMatrix
        IntegerCharacterMatrix matrix = value.value();
        // Create a new BEAST Alignment object that will hold the copy number data
        Alignment alignment = new Alignment();
        // Initialize a list to store the sequence data for each taxon
        List<Sequence> sequences = new ArrayList<>();

        // Calculate maximum copy number
        int maxCopyNumber = 0;
        for (int i = 0; i < matrix.getTaxa().ntaxa(); i++) {
            String taxonName = String.valueOf(matrix.getTaxa().getTaxaNames()[i]);
            for (int j = 0; j < matrix.nchar(); j++) {
                int state = matrix.getState(taxonName, j);
                if (state > maxCopyNumber) {
                    maxCopyNumber = state;
                }
            }
        }
        int totalCount = maxCopyNumber + 1;

        // Iterate through each taxon in the matrix
        for (int i = 0; i < matrix.getTaxa().ntaxa(); i++) {
            String taxonName = String.valueOf(matrix.getTaxa().getTaxaNames()[i]);
            context.addTaxon(taxonName);

            StringBuilder sequenceData = new StringBuilder();
            for (int j = 0; j < matrix.nchar(); j++) {
                // Add comma between values (but not before the first value)
                if (j > 0) {
                    sequenceData.append(",");
                }
                // Get and append the integer state
                int state = matrix.getState(taxonName, j);
                sequenceData.append(state);
            }
            // Create a new BEAST sequence
            Sequence sequence = new Sequence();
            sequence.setInputValue("taxon", taxonName);
            sequence.setInputValue("value", sequenceData.toString());
            sequence.setInputValue("totalcount", totalCount);
            sequence.initAndValidate();
            // Add the sequence to 'sequence' list
            sequences.add(sequence);
        }
        // Add all sequences to the alignment
        alignment.setInputValue("sequence", sequences);
        // Set the data type to integer
        alignment.setInputValue("dataType", "integer");

        // Initialize the alignment
        alignment.initAndValidate();
        // Return the completed alignment
        return alignment;
    }


    @Override
    public Class getValueClass() {
        return IntegerCharacterMatrix.class;
    }

    @Override
    public Class<Alignment> getBEASTClass() {
        return Alignment.class;
    }
}
