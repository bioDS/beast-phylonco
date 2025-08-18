package phylonco.lphy.evolution.copynumbermodel;

import java.io.BufferedReader;

// this class should read in copy number matrix data from txt files

@IOFunction(
        role = IOFunction.Role.dataInput,
        extensions = { ".txt"},
        fileArgument = ReaderConst.FILE
)
public class ReadCopyProfile extends DeterministicFunction<IntegerCharacterMatrix> {

	public ReadFasta(@ParameterInfo(name = ReaderConst.FILE, description = "the name of txt file of copy number profiles.") Value<String> filePath) {

        if (filePath == null) {
        	throw new IllegalArgumentException("The file name can't be null!");
        }

        setParam(ReaderConst.FILE, filePath);

    }


    public Value<IntegerCharacterMatrix> apply() { 
    	String filePath = ((Value<String>) getParams().get(ReaderConst.FILE)).value();

    	Path readPath = UserDir.getUserPath(filePath);

    	// note: bufferedreader may throw an exception
    	// hint: add a try and catch block here to fail elegantly
    	BufferedReader reader = new BufferedReader(new FileReader(readPath));

    	IntegerCharacterMatrix data = getData(reader);

    	return new Value<>(null, data, this);

    }

	private IntegerCharacterMatrix getData(BufferedReader reader) { 
		// TODO: implement this method to parse data from BufferedReader input

		Taxa taxa = null; // read taxa names and construct Taxa objects

		int bins = 100; // get numnber of bins from txt input file

		IntegerCharacterMatrix matrix = new IntegerCharacterMatrix(taxa, bins);

		// hint: you'll need a for loop here
		int taxonIndex = 0;
		int position = 0;
		int state = 0;

		// set state for each bin and taxa
		matrix.setState(taxonIndex, position, state);

		return matrix;
	}


}
