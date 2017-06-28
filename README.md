ALE
[![Build Status](https://travis-ci.org/ssolo/ALE.svg?branch=master)](https://travis-ci.org/ssolo/ALE)
===
## ALE: Reconstruction of reconciled gene trees and estimation of DTL rates from a sample of gene trees and a rooted species tree.

#### Installation instructions are provided at the bottom

Amalgamated likelihood estimation (ALE) is a probabilistic approach to exhaustively explore all reconciled gene trees that can be amalgamated as a combination of clades observed in a sample of gene trees. We implement the ALE approach in the context of a reconciliation model (cf. http://arxiv.org/abs/1211.4606 ), which allows for the duplication, transfer and loss of genes. We use ALE to efficiently approximate the sum of the joint likelihood over amalgamations and to find the reconciled gene tree that maximizes the joint likelihood among all such trees.   

The repository currently implements two objects and several programs:

1. The fundamental approx_posterior object that can be used to observe and store split frequencies and provides interfaces for performing basic Conditional Clade Probability (CCP, cf. Hohna and Drummond Systbiol. 2012) related, as well as caluclations related to more general Amalgamated Likelihood Estimation (ALE, cf. "Efficient Exploration of the Space of Reconciled Gene Trees", Syst Biol. 2013 Nov;62(6):901-12. doi: 10.1093/sysbio/syt054.).  

2. The exODT_model object that implements ALE in the context of our recently published reconciliation model (cf. http://arxiv.org/abs/1211.4606 ).

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

#### Bayesian MCMC sampling with ALEmcmc_undated
DTL rates as well as the origination at the root rate are sampled. A prior value for these rates can be provided:
```sh
ALEmcmc_undated species_tree.newick geneFamily.ale sample=number_of_samples separators=gene_name_separator O_R=OriginationAtRootPrior delta=DuplicationRatePrior tau=TransferRatePrior lambda=LossRatePrior sampling_rate=sampling_rate beta=weight_of_sequence_evidence
```
All the parameters correspond to those used by ALEml_undated, except for :

sampling_rate=sampling_rate: during MCMC, sample a tree every sampling_rate iterations.

For more documentation please consult our manuscript "Efficient Exploration of the Space of Reconciled Gene Trees" (Syst Biol. 2013 Nov;62(6):901-12. doi: 10.1093/sysbio/syt054.) together with the code and comments for the ALEobserve, ALEml_undated, ALEmcmc_undated, ALEml, ALEmcmc programs that deploy these objects.

Example data and an associated README file are available in the example_data subdirectory.   


## Installation instructions

There are two ways to use the ALE suite of programs. The simplest approach is using Docker, the other uses cmake and requires installing several dependencies.

### Simple Docker installation

To use ALE using Docker, you need to install Docker first.
* Under Mac OS X, this is done with a DMG file that can be downloaded from [there](https://docs.docker.com/docker-for-mac/install/#install-and-run-docker-for-mac).
* Under Windows, you can download the executable file from [there](https://docs.docker.com/docker-for-windows/install/)
* Under Linux, you need to choose your particular flavor on the left of [this page](https://docs.docker.com/engine/installation/) and then follow the instructions.

Once Docker has been installed, you can use ALE as follows. As described above, you need to run ALEobserve then ALEml or ALEmcmc_undated to get reconciled gene trees, and this is explained below.

Let's assume you have a file named "geneFamily.treelist" containing a distribution of gene trees, and a species tree named species_tree.newick.
We will launch the program from the folder containing those data, which means that we will launch the command from the folder $PWD. $PWD is the usual environmental variable that UNIX-type operating systems use to point to the Present Working Directory.
If you want to use the example data to test the software, this means your $PWD will be /absolute_path/ALE/example_data, your species tree will be S.tree, and your gene tree file will be either HBG745965_real.1.treelist.txt or HBG745965_real.2.treelist.txt.

The way the docker image is as follows: you can run all the programs of the ALE suite, and only have to precede the command you would use to launch them with "docker run -v $PWD:$PWD alesuite ".

For instance, we provide below what the typical pipeline would look like.


- The first command to use is:
```sh
docker run -v $PWD:$PWD alesuite ALEobserve $PWD/geneFamily.treelist
```

- The second command to use could be:
```sh
docker run -v $PWD:$PWD alesuite ALEml_undated $PWD/species_tree.newick $PWD/geneFamily.treelist.ale
```

or

```sh
docker run -v $PWD:$PWD alesuite ALEmcmc_undated $PWD/species_tree.newick $PWD/geneFamily.treelist.ale
```

for instance.

Running ALEobserve then one of those two commands will produce the output files describing scenarios of gene family evolution including events of gene transfer, duplication and loss.


### Expert installation using cmake

WARNING: We have recently transitioned to cmake and installation may be buggy. If you encounter any problem please let us know.  

We have successfully installed ALE on Unix systems (Linux and Mac OS X systems).

Installation using cmake requires:
 - cmake
 - a C++ compiler (e.g. g++ or clang)
 - the Bio++ libraries bpp-core, bpp-seq, and bpp-phyl
 - the Boost C++ libraries

We advise that you create a folder named build in the ALE folder:

```sh
mkdir build
cd build
```

Then run cmake:

```sh
cmake ..
make
```

Using make with the option "-j 4" paralellizes compilation across 4 processes and speeds it up. Those commands will produce executable files in the folder ALE/build/bin.


All code, data and documentation by Szollosi GJ et al.; ssolo@elte.hu; CC BY-SA 3.0;
