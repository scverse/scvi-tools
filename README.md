# scVI

### Title
A deep generative model for gene expression profiles
from single-cell RNA sequencing

### Authors
Romain Lopez, Jeffrey Regier, Michael Cole, Michael Jordan and Nir Yosef <br />
Department of Electrical Engineering and Computer Sciences <br />
Department of Physics <br />
University of California, Berkeley <br />

### Abstract
Transcriptome profiles of single cells reflect true and often unexplored biological diversity, but are also affected by noise of biological and technical nature. This raises the need to explicitly model the resulting uncertainty and take it into account in any downstream analysis. Here, we introduce scVI (single cell variational inference) - a  framework for probabilistic interpretation of gene expression in single cells. Our model uses a variational autoencoder to approximate the parameters that govern the distribution of expression values of each gene in every cell, using a non linear mapping of the observations into a low-dimensional latent space. By doing so, scVI pools information between similar cells and genes while taking confounding factors such as batch, varying dropout rates and library size into account. We conduct a comprehensive comparative analysis to existing methods for distributional modeling and dimensionality reduction, all of which rely on generalized linear models. We first show that the scVI inference procedure scales to over one million cells, whereas competing algorithms do not. Next, we show that scVI fits unseen data more closely and can impute missing data more accurately, both indicative of a better generalization capacity. We then extend scVI to conduct basic analytical tasks such as batch removal, visualization, clustering and differential expression and demonstrate its accuracy in comparison to the state of the art tools in each task. scVI is publicly available, and can be readily used as a principled and inclusive solution for multiple tasks of single cell RNA-seq data analysis.

## Prerequisites
+ Python 2.7
+ TensorFlow 1.2.0

## Codes
In this github repo, we will provide you with code to:
+ minimal code to run the model
+ a complete walk through the model architecture
+ reproduce results of the scVI paper

## Additional comments
Some effort might be needed to have a functional installation of the GPU version of TensorFlow. A great tutorial is given here: https://www.tensorflow.org/install/install_linux#gpu_support

### Contact
romain_lopez@berkeley.edu
