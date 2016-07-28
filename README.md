ALE
===
ALE: Reconstruction of reconciled gene trees and estimation of DTL rates from a sample of gene trees and a rooted species tree.


WARNING: We have recently transitioned to cmake and installation may be buggy. If you encounter any problem please let us know.   

Amalgamated likelihood estimation (ALE) is a probabilistic approach to exhaustively explore all reconciled gene trees that can be amalgamated as a combination of clades observed in a sample of gene trees. We implement the ALE approach in the context of a reconciliation model (cf. http://arxiv.org/abs/1211.4606 ), which allows for the duplication, transfer and loss of genes. We use ALE to efficiently approximate the sum of the joint likelihood over amalgamations and to find the reconciled gene tree that maximizes the joint likelihood among all such trees.   

The repository currently implements two objects and several programs:

1. The fundamental approx_posterior object that can be used to observe and store split frequencies and provides interfaces for performing basic Conditional Clade Probability (CCP, cf. Hohna and Drummond Systbiol. 2012) related, as well as caluclations related to more general Amalgamated Likelihood Estimation (ALE, cf. "Efficient Exploration of the Space of Reconciled Gene Trees", Syst Biol. 2013 Nov;62(6):901-12. doi: 10.1093/sysbio/syt054.).  

2. The exODT_model object that implements ALE in the context of our recently published reconciliation model (cf. http://arxiv.org/abs/1211.4606 ).

## Using ALE
ALE reconciles a sample of gene trees with a species tree. The sample of gene trees is typically obtained from a Bayesian MCMC program (e.g. PhyloBayes or mrBayes or RevBayes...) or could be obtained using bootstrap replicates. The species tree needs to be rooted. Both inputs need to be provided as files containing Newick strings. ALE will use the sample of gene trees to piece together (amalgamate) reconciled gene trees while estimating their probabilities. Here the probability of a gene tree depends on both sequence information (provided by the gene tree frequencies in the initial sample) and rates of D, T, L events, which can be optimized or sampled depending on the program chosen. 2 models have been implemented for reconciling gene trees with a species tree. A dated model, which assumes that the nodes of the species tree are ordered relative to each other, and an undated model which only requires that the species tree is rooted.

Two steps are required to obtained reconciled amalgamated gene trees.
1. First, compact the sample of gene trees into an ALE object. This is done using the program ALEobserve.
```sh
ALEobserve geneFamily.treelist 1000  
```
Where geneFamily.treelist is the file containing the sample of gene trees, and 1000 is a number of gene trees that will be discarded. The result of the command above is to create a geneFamily.ale file that can be used to infer reconciled amalgamated gene trees during step 2.

2. Second, compute the reconciled amalgamated gene trees. The DTL rates can be provided, or estimated. For the undated model, there are 2 ways to estimate these rates: Maximum Likelihood estimation with ALEml_undated, or Bayesian MCMC sampling with ALEmcmc_undated.

### Maximum likelihood estimation with ALEml_undated
DTL rates can either be estimated by providing no starting value or fixed to a starting value. Using those rates, the program then outputs a certain number of sampled reconciled amalgamated gene trees.
```sh
ALEml_undated species_tree.newick geneFamily.ale sample=number_of_samples separators=gene_name_separator O_R=OriginationAtRoot delta=DuplicationRate tau=TransferRate lambda=LossRate beta=weight_of_sequence_evidence
```
species_tree.newick : contains the Newick description of the rooted species tree.
geneFamily.ale : contains the ALE built during step 1.
sample=number_of_samples : how many reconciled gene trees should be output.
separators=gene_name_separator : separator character used to separate species names from gene names (default to "_").
O_R=OriginationAtRoot : rate of origination at root.
delta=DuplicationRate : duplication rate.
tau=TransferRate : transfer rate.
lambda=LossRate : loss rate.
beta=weight_of_sequence_evidence : how much sequence data will be trusted. Values higher than 1 mean the gene trees coming from the sequences alone is more trusted than by default; values lower than 1 mean that the gene trees seem not very trustworthy. Defaults to 1.

### Bayesian MCMC sampling with ALEmcmc_undated
DTL rates as well as the origination at the root rate are sampled. A prior value for these rates can be provided:
```sh
ALEmcmc_undated species_tree.newick gene_tree_sample.ale sample=number_of_samples separators=gene_name_separator O_R=OriginationAtRootPrior delta=DuplicationRatePrior tau=TransferRatePrior lambda=LossRatePrior sampling_rate=sampling_rate beta=weight_of_sequence_evidence
```
All the parameters correspond to those used by ALEml_undated, except for :
sampling_rate=sampling_rate: during MCMC, sample a tree every sampling_rate iterations.

For more documentation please consult our manuscript "Efficient Exploration of the Space of Reconciled Gene Trees" (Syst Biol. 2013 Nov;62(6):901-12. doi: 10.1093/sysbio/syt054.) together with the code and comments for the ALEobserve, ALEml_undated, ALEmcmc_undated, ALEml, ALEmcmc programs that deploy these objects.

Example data and an associated README file are available in the example_data subdirectory.   

All code, data and documentation by Szollosi GJ et al.; ssolo@elte.hu; CC BY-SA 3.0;
