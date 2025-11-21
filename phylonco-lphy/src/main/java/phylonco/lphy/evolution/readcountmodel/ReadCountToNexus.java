package phylonco.lphy.evolution.readcountmodel;

import lphy.core.model.DeterministicFunction;
import lphy.core.model.Value;
import lphy.core.model.annotation.GeneratorInfo;
import lphy.core.model.annotation.ParameterInfo;

public class ReadCountToNexus extends DeterministicFunction<ReadCountNexus> {
    public ReadCountToNexus(@ParameterInfo(name = "rc", description = "read count data that converts to nexus file") Value<ReadCountData> rc) {
        if (rc == null) {
            throw new IllegalArgumentException("Read count data should not be null");
        }
        setParam("rc", rc);
    }

    @GeneratorInfo(name = "toNexus", examples = {"mpileupToReadCount.lphy"},
            description = "covert the read count data to nexus file")
    @Override
    public Value<ReadCountNexus> apply() {
        ReadCountData rc = getReadCountData().value();
        ReadCountNexus nexus = new ReadCountNexus(rc);

        return new Value<>("", nexus, this);
    }

    public Value<ReadCountData> getReadCountData() {
        return getParams().get("rc");
    }
}
