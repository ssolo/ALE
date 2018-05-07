ALE
[![Build Status](https://travis-ci.org/ssolo/ALE.svg?branch=master)](https://travis-ci.org/ssolo/ALE)
===
## ALE: Reconstruction of reconciled gene trees and estimation of DTL rates from a sample of gene trees and a rooted species tree.

#### For installation instructions, please refer to the [INSTALL file](INSTALL.md).


Amalgamated likelihood estimation (ALE) is a probabilistic approach to exhaustively explore all reconciled gene trees that can be amalgamated as a combination of clades observed in a sample of gene trees. We implement the ALE approach in the context of a reconciliation model (for ALE dated, cf. [our Syst. Biol. paper](http://arxiv.org/abs/1211.4606) and for ALE_undated, cf. [our Phil. Trans. B paper](http://rstb.royalsocietypublishing.org/content/370/1678/20140335.long)), which allows for the duplication, transfer and loss of genes. We use ALE to efficiently approximate the sum of the joint likelihood over amalgamations and to find the reconciled gene tree that maximizes the joint likelihood among all such trees.   

The repository currently implements two objects and several programs:

1. The fundamental approx_posterior object that can be used to observe and store split frequencies and provides interfaces for performing basic Conditional Clade Probability (CCP, cf. [Hohna and Drummond Syst. Biol. 2012](https://academic.oup.com/sysbio/article/61/1/1/1676649/Guided-Tree-Topology-Proposals-for-Bayesian)) related, as well as calculations related to more general Amalgamated Likelihood Estimation (ALE, cf. ["Efficient Exploration of the Space of Reconciled Gene Trees", Syst Biol. 2013 Nov;62(6):901-12. doi: 10.1093/sysbio/syt054.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3797637/)).  

2. The exODT_model object that implements ALE in the context of our recently published reconciliation model (for ALE dated, cf. [our Syst. Biol. paper](http://arxiv.org/abs/1211.4606) and for ALE_undated, cf. [our Phil. Trans. B paper](http://rstb.royalsocietypublishing.org/content/370/1678/20140335.long)).

## Using ALE
ALE reconciles a sample of gene trees with a species tree. The sample of gene trees is typically obtained from a Bayesian MCMC program (e.g. PhyloBayes or mrBayes or RevBayes...) or may be obtained using bootstrap replicates (although this would be less correct). The species tree needs to be rooted. Both inputs need to be provided as files containing Newick strings. ALE will use the sample of gene trees to piece together (amalgamate) reconciled gene trees while estimating their probabilities. Here the probability of a gene tree depends on both sequence information (provided by the gene tree frequencies in the initial sample) and rates of D, T, L events, which can be optimized or sampled depending on the program chosen. 2 models have been implemented for reconciling gene trees with a species tree. A dated model, which assumes that the nodes of the species tree are ordered relative to each other, and an undated model which only requires that the species tree is rooted.

Two steps are required to obtained reconciled amalgamated gene trees: first the preparation of the data by constructing an ALE object, then the actual inference of the reconciled amalgamated gene trees.

### Construction of the ALE object
The sample of gene trees is compacted into an ALE object. This is done using the program ALEobserve.
```sh
ALEobserve geneFamily.treelist 1000  
```
Where geneFamily.treelist is the file containing the sample of gene trees, and 1000 is a number of gene trees that will be discarded. The result of the command above is to create a geneFamily.ale file that can be used to infer reconciled amalgamated gene trees during step 2.


### Inference of the reconciled gene trees

The second step infers the reconciled amalgamated gene trees. The DTL rates can be provided, or estimated. For the undated model, there are 2 ways to estimate these rates: Maximum Likelihood estimation with ALEml_undated, or Bayesian MCMC sampling with ALEmcmc_undated.

#### Maximum likelihood estimation with ALEml_undated
DTL rates can either be estimated by providing no starting value or fixed to a starting value. Using those rates, the program then outputs a certain number of sampled reconciled amalgamated gene trees.
```sh
ALEml_undated species_tree.newick geneFamily.ale sample=number_of_samples separators=gene_name_separator O_R=OriginationAtRoot delta=DuplicationRate tau=TransferRate lambda=LossRate beta=weight_of_sequence_evidence fraction_missing=fraction_missing.txt
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

fraction_missing=fraction_missing.txt : file containing the expected fraction of missing genes per species, in the format "species_name:fraction", where species_name should match the species name in the species tree, and fraction should be a floating number between 0 and 1. There should be as many lines as there are species.

#### Bayesian MCMC sampling with ALEmcmc_undated
DTL rates as well as the origination at the root rate are sampled. A prior value for these rates can be provided:
```sh
ALEmcmc_undated species_tree.newick geneFamily.ale sample=number_of_samples separators=gene_name_separator O_R=OriginationAtRootPrior delta=DuplicationRatePrior tau=TransferRatePrior lambda=LossRatePrior sampling_rate=sampling_rate beta=weight_of_sequence_evidence
```
All the parameters correspond to those used by ALEml_undated, except for :

sampling_rate=sampling_rate: during MCMC, sample a tree every sampling_rate iterations.

For more documentation please consult our manuscript "Efficient Exploration of the Space of Reconciled Gene Trees" (Syst Biol. 2013 Nov;62(6):901-12. doi: 10.1093/sysbio/syt054.) together with the code and comments for the ALEobserve, ALEml_undated, ALEmcmc_undated, ALEml, ALEmcmc programs that deploy these objects.

Example data and an associated README file are available in the example_data subdirectory.   


All code, data and documentation by Szollosi GJ et al.; ssolo@elte.hu;

    This file is part of ALE.

    ALE is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    ALEr is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Foobar.  If not, see <http://www.gnu.org/licenses/>.

