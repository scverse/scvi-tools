"""Perform the sanity checks once deconvolution type 1 is performed."""
# flake8: noqa
# %%
import scanpy as sc
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import random
import mygene
from sklearn.model_selection import train_test_split
from scipy.stats import pearsonr

from sanity_checks_utils import (
    GROUPS,
    read_almudena_signature,
    map_hgnc_to_ensg,
    perform_nnls,
    compute_correlations,
)

# %%
# import data
adata = sc.read("/home/owkin/data/cross-tissue/omics/raw/local.h5ad")
signature_laughney = pd.read_csv(
    "/home/owkin/data/laughney_signature.csv", index_col=0
).drop(["Endothelial", "Malignant", "Stroma", "Epithelial"], axis=1)
signature_almudena = read_almudena_signature(
    "/home/owkin/project/Almudena/Output/Crosstiss_Immune/CTI.txt"
)
# primary cell type categories
groups = GROUPS["primary_groups"]


# %%
# create cell types
group_correspondence = {}
for k, v in groups.items():
    for cell_type in v:
        group_correspondence[cell_type] = k
adata.obs["cell_types_grouped"] = [
    group_correspondence[cell_type] for cell_type in adata.obs.Manually_curated_celltype
]

# remove some cell types: you need more than 15GB memory to run that
index_to_keep = adata.obs.loc[adata.obs["cell_types_grouped"] != "To remove"].index
adata = adata[index_to_keep]


# %%
# build signature on train set and apply deconvo on the test set
cell_types_train, cell_types_test = train_test_split(
    adata.obs_names,
    test_size=0.5,
    stratify=adata.obs.cell_types_grouped,
    random_state=42,
)
adata = adata[cell_types_test, :]


# %%

signature_choice = "almudena"  # almudena or laughney

if signature_choice == "laughney":
    # map the HGNC notation to ENSG if the signature matrix uses HGNC notation
    mg = mygene.MyGeneInfo()
    genes = mg.querymany(
        signature_laughney.index,
        scopes="symbol",
        fields=["ensembl"],
        species="human",
        verbose=False,
        as_dataframe=True,
    )
    ensg_names = map_hgnc_to_ensg(genes, adata)
    signature = signature_laughney.copy()
    signature.index = ensg_names
elif signature_choice == "almudena":
    signature = signature_almudena.copy()

# intersection between all genes and marker genes
intersection = list(set(adata.var_names).intersection(signature.index))
signature = signature.loc[intersection]


# %%
# Sanity check 1: Pseudobulks for each cell type → the closest the fraction is to 1 the best

# The following is the equivalent of averaged_data = adata.groupby("cell_types_grouped").X.mean(axis = 0)
# but there is no "groupby" function for an anndata sparse matrix
grouped = adata.obs.groupby("cell_types_grouped")
averaged_data, group = [], []
for group_key, group_indices in grouped.groups.items():
    averaged_data.append(adata[group_indices].X.mean(axis=0).tolist()[0])
    group.append(group_key)

averaged_data = pd.DataFrame(averaged_data, index=group, columns=adata.var_names)

# intersection between all genes and marker genes
averaged_data = averaged_data[intersection]

# 5 different deconvolutions with n_marker_genes observations each
deconv_results = perform_nnls(signature, averaged_data)

# melt the matrix for seaborn
deconv_results_melted = pd.melt(
    deconv_results.T.reset_index(),
    id_vars="index",
    var_name="Cell type",
    value_name="Fraction",
).rename({"index": "Cell type predicted"}, axis=1)

# if one deconv method only
sns.set_style("whitegrid")
sns.stripplot(
    data=deconv_results_melted, x="Cell type", y="Fraction", hue="Cell type predicted"
)
plt.show()

# if comparison of several deconv methods
deconv_results_melted_methods = deconv_results_melted.loc[
    deconv_results_melted["Cell type predicted"] == deconv_results_melted["Cell type"]
].copy()
deconv_results_melted_methods["Method"] = "nnls"  # add all the deconv methods
sns.set_style("whitegrid")
sns.stripplot(
    data=deconv_results_melted_methods, x="Cell type", y="Fraction", hue="Method"
)
plt.show()


# %%
# Sanity check 2: Pseudobulks from predefined fractions → correlation between the
# fractions of selected cells (“ground truth”)  and the estimated fractions

n_sample = 100
random.seed(42)
averaged_data = []
ground_truth_fractions = []
for i in range(n_sample):
    cell_sample = random.sample(list(adata.obs_names), 1000)
    adata_sample = adata[cell_sample, :]
    ground_truth_frac = adata_sample.obs.cell_types_grouped.value_counts() / 1000
    ground_truth_fractions.append(ground_truth_frac)
    averaged_data.append(adata_sample.X.mean(axis=0).tolist()[0])

averaged_data = pd.DataFrame(
    averaged_data, index=range(n_sample), columns=adata.var_names
)
ground_truth_fractions = pd.DataFrame(ground_truth_fractions, index=range(n_sample))
ground_truth_fractions = ground_truth_fractions.fillna(
    0
)  # the Nan are cells not sampled

# intersection between all genes and marker genes
averaged_data = averaged_data[intersection]

# 100 different deconvolutions with n_marker_genes observations each
deconv_results = perform_nnls(signature, averaged_data)

# compute correlations
correlations = compute_correlations(deconv_results, ground_truth_fractions)

# boxplot of the correlations for the method
sns.set_style("whitegrid")
boxplot = sns.boxplot(correlations, y="correlations", x="Method")
medians = correlations.groupby(["Method"])["correlations"].median().round(4)
vertical_offset = (
    correlations["correlations"].median() * 0.0005
)  # for non-overlapping display
for xtick in boxplot.get_xticks():  # show the median value
    boxplot.text(
        xtick,
        medians[xtick] + vertical_offset,
        medians[xtick],
        horizontalalignment="center",
        size="x-small",
        color="w",
        weight="semibold",
    )
plt.show()


# %%
# Sanity check 3: Pseudobulks from random selection of cells (e.g. Dirichlet
# distribution) → correlation between the fractions of selected cells (“ground truth”)
# and the estimated fractions

n_sample = 100

# compute dirichlet posteriors to sample cells
# nb: dirichlet is conjugate to the multinomial distribution, thus giving an easy
# posterior calculation
np.random.seed(42)
prior_alphas = np.ones(
    len(signature.columns)
)  # non-informative prior: if prior belief, should be incorporated here
likelihood_alphas = adata.obs.cell_types_grouped.value_counts() / len(
    adata.obs
)  # multinomial likelihood (TO CHECK: maths behind average or sum for the likelihood)
alpha_posterior = (
    prior_alphas + likelihood_alphas
)  # (TO CHECK: maths behind this simple addition)
posterior_dirichlet = np.random.dirichlet(alpha_posterior, n_sample)
posterior_dirichlet = np.round(posterior_dirichlet * 500)
posterior_dirichlet = posterior_dirichlet.astype(np.int64)  # number of cells to sample

ground_truth_fractions = posterior_dirichlet / posterior_dirichlet.sum(
    axis=1, keepdims=True
)
ground_truth_fractions = pd.DataFrame(
    ground_truth_fractions, columns=alpha_posterior.index
)

random.seed(42)
averaged_data = []
for i in range(n_sample):
    sample_data = []
    for j, cell_type in enumerate(likelihood_alphas.index):
        # is there a way to optimise by not doing this loop?
        cell_sample = random.sample(
            list(adata.obs.loc[adata.obs.cell_types_grouped == cell_type].index),
            posterior_dirichlet[i][j],
        )
        adata_sample = adata[cell_sample, :]
        sample_data.append(adata_sample.X.mean(axis=0).tolist()[0])
    averaged_data.append(np.array(sample_data).mean(axis=0))

averaged_data = pd.DataFrame(
    averaged_data, index=range(n_sample), columns=adata.var_names
)

# intersection between all genes and marker genes
averaged_data = averaged_data[intersection]

# 100 different deconvolutions with n_marker_genes observations each
deconv_results = perform_nnls(signature, averaged_data)

# compute correlations
correlations = compute_correlations(deconv_results, ground_truth_fractions)

# boxplot of the correlations for the method
sns.set_style("whitegrid")
boxplot = sns.boxplot(correlations, y="correlations", x="Method")
medians = correlations.groupby(["Method"])["correlations"].median().round(4)
vertical_offset = (
    correlations["correlations"].median() * 0.0005
)  # for non-overlapping display
for xtick in boxplot.get_xticks():  # show the median value
    boxplot.text(
        xtick,
        medians[xtick] + vertical_offset,
        medians[xtick],
        horizontalalignment="center",
        size="x-small",
        color="w",
        weight="semibold",
    )
plt.show()

# %%
