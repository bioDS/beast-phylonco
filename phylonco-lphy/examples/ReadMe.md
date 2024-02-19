### Newly Added File Explanation
Three LPhy scripts newly added need explanations:

1. ErrModel_Haploid.lphy
2. true_alignment.lphy
3. test_ambiguities.lphy

In the experiment, we consider part of the human reference genome as the root sequence to 
generate coalescent trees with GT16 substitution model and obtain the true alignment. We assume
the genotype is homozygous, so that we write a new Generative Distribution function called 
"Homozygous" to convert the input alignment into homozygous alignment. Meanwhile, if there are
ambiguities or gap in the input sequence, the function would randomly pick a possible state on 
that site for generating the alignment. To sequence the alignments, haploid fasta files are 
needed. Therefore, a new deterministic function "haploid" was written to split the diploid 
alignment into two haploid alignments. 

The first two scripts are used for simulating sequences in the experiments. "true_alignment.lphy"
simulates the true alignment while "ErrModel_Haploid.lphy" simulates the alignment corrected by
the GT16 error model. "test_ambiguities.lphy" is used to test the ambiguities and gap processing.

