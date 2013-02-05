ALE
===

Amalgamated likelihood estimation (ALE) is a probabilistic approach to exhaustively explore all reconciled gene trees that can be amalgamated as a combination of clades observed in a sample of gene trees. We implement the ALE approach in the context of a reconciliation model (cf. http://arxiv.org/abs/1211.4606 ), which allows for the duplication, transfer and loss of genes. We use ALE to efficiently approximate the sum of the joint likelihood over amalgamations and to find the reconciled gene tree that maximizes the joint likelihood among all such trees.   

The repository currently impliemnets two objects:

1. The fundamental approx_posterior object that can be used to observe and store split frequencies and provides interfaces for performing basic Conditional Clade Probability (CCP, cf. Hohna and Drummond Systbiol. 2012) related, as well as more general Amalgamated Likelihood Estimation (ALE, cf. "Efficient Exploration of the Space of Reconciled Gene Trees", under review, soon available from the arxiv) related calculations.  

2. The exODT_model object that impliments ALE in the context of our recently published reconcilation model (cf. http://arxiv.org/abs/1211.4606 ).


In lieu of proper documentation please consult our manuscript "Efficient Exploration of the Space of Reconciled Gene Trees" (under review, soon available from the arxiv) together with the code and comments for the ALEobserve, ALEml, ALEsample programs that deploy these objects. 

Example data and an associated README file are available in the example_data subdirectory.   

Precompiled binaries for OSX and LINUX (including OpenMP versions) are available in the binaries subdirectory. 
