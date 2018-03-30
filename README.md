# single-cell Variational Inference (scVI)

+ NIPS MLCB Submission: https://arxiv.org/abs/1709.02082
+ bioRxiv preprint: https://www.biorxiv.org/content/early/2018/03/30/292037
+ implementation for the preprint: https://github.com/romain-lopez/scVI-reproducibility

### Title
Bayesian Inference for a Generative Model of Transcriptome Profiles from Single-cell RNA Sequencing

### Authors
Romain Lopez, Jeffrey Regier, Michael Cole, Michael Jordan and Nir Yosef <br />
Department of Electrical Engineering and Computer Sciences <br />
Department of Physics <br />
University of California, Berkeley <br />

### Abstract
Transcriptome profiles of individual cells reflect true and often unexplored biological diversity, but are also affected by noise of biological and technical nature. This raises the need to explicitly model the resulting uncertainty and take it into account in any downstream analysis, such as dimensionality reduction, clustering, and differential expression. Here, we introduce Single-cell Variational Inference (scVI), a scalable framework for probabilistic representation and analysis of gene expression in single cells. Our model uses variational inference and stochastic optimization of deep neural networks to approximate the parameters that govern the distribution of expression values of each gene in every cell, using a non-linear mapping between the observations and a low-dimensional latent space.


By doing so, scVI pools information between similar cells or genes while taking nuisance factors of variation such as batch effects and limited sensitivity into account. To evaluate scVI, we conducted a comprehensive comparative analysis to existing methods for distributional modeling and dimensionality reduction, all of which rely on generalized linear models. We first show that scVI scales to over one million cells, whereas competing algorithms can process at most tens of thousands of cells. Next, we show that scVI fits unseen data more closely and can impute missing data more accurately, both indicative of a better generalization capacity. We then utilize scVI to conduct a set of fundemental analysis tasks -- including batch correction, visualization, clustering and differential expression -- and demonstrate its accuracy in comparison to the state-of-the-art tools in each task. scVI is publicly available, and can be readily used as a principled and inclusive solution for multiple tasks of single-cell RNA sequencing data analysis.

## Prerequisites
+ Python 2.7
+ TensorFlow 1.2.0

## Codes
In this github repo, we will provide you with minimal code to run the model on a public dataset

## Additional comments
Some effort might be needed to have a functional installation of the GPU version of TensorFlow. A great tutorial is given here: https://www.tensorflow.org/install/install_linux#gpu_support

### Contact
romain_lopez@berkeley.edu
