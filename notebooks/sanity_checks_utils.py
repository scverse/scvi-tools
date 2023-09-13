"""Different python functions useful for sanity checks in deconvolution."""

import pandas as pd
import numpy as np
import scipy.stats
from sklearn.linear_model import LinearRegression

GROUPS = {
    "primary_groups": {
        "B": [
            "ABCs",
            "GC_B (I)",
            "GC_B (II)",
            "Memory B cells",
            "Naive B cells",
            "Plasma cells",
            "Plasmablasts",
            "Pre-B",
            "Pro-B",
        ],
        "MonoMacro": [
            "Alveolar macrophages",
            "Classical monocytes",
            "Erythrophagocytic macrophages",
            "Intermediate macrophages",
            "Nonclassical monocytes",
        ],
        "TNK": [
            "Cycling T&NK",
            "MAIT",
            "NK_CD16+",
            "NK_CD56bright_CD16-",
            "T_CD4/CD8",
            "Teffector/EM_CD4",
            "Tem/emra_CD8",
            "Tfh",
            "Tgd_CRTAM+",
            "Tnaive/CM_CD4",
            "Tnaive/CM_CD4_activated",
            "Tnaive/CM_CD8",
            "Tregs",
            "Trm/em_CD8",
            "Trm_Tgd",
            "Trm_Th1/Th17",
            "Trm_gut_CD8",
            "ILC3",
        ],
        "DC": ["DC1", "DC2", "migDC", "pDC"],
        "Mast": ["Mast cells"],
        "To remove": [
            "Erythroid",
            "Megakaryocytes",
            "Progenitor",
            "Cycling",
            "T/B doublets",
            "MNP/B doublets",
            "MNP/T doublets",
            "Intestinal macrophages",
        ],
    },
    "precise_groups": {
        "B": [
            "ABCs",
            "GC_B (I)",
            "GC_B (II)",
            "Memory B cells",
            "Naive B cells",
            "Pre-B",
            "Pro-B"
        ],
        "Plasma": ["Plasma cells", "Plasmablasts"],
        "Mono": ["Classical monocytes", "Nonclassical monocytes"],
        "CD8 T": ["Tem/emra_CD8", "Tnaive/CM_CD8", "Trm/em_CD8", "Trm_gut_CD8"],
        "CD4 T":["Teffector/EM_CD4", "Tfh", "Tnaive/CM_CD4", "Tnaive/CM_CD4_activated", "Tregs", 
                 "Trm_Th1/Th17"],
        "T": ["MAIT","T_CD4/CD8","Tgd_CRTAM+","Trm_Tgd"],
        "NK": ["NK_CD16+", "NK_CD56bright_CD16-"],
        "DC": ["DC1", "DC2", "migDC", "pDC"],
        "Mast": ["Mast cells"],
        "Red_blood": ["Erythroid"],
        "Bone_marrow": ["Megakaryocytes"],
        "Non-differentiated": ["Progenitor"],
        "To remove":  ["Cycling", "T/B doublets", "Cycling T&NK",
            "MNP/B doublets", "MNP/T doublets","Alveolar macrophages",
            "Erythrophagocytic macrophages",
            "Intermediate macrophages",
            "Intestinal macrophages", "ILC3"],
    },
}


def read_almudena_signature(path):
    """Read Almudena's signature matrix. Requires this function because it's a txt file
    delimited with various delimiters.
    """
    signature_almudena = []
    with open(path, "r") as file:
        for line in file:
            temp = []
            for elem in line.split("\t"):
                try:
                    temp.append(float(elem))
                except:
                    elem = elem.replace('"', "")
                    elem = elem.replace("\n", "")
                    temp.append(elem)
            signature_almudena.append(temp)

    signature_almudena = pd.DataFrame(signature_almudena).set_index(0)
    signature_almudena.columns = signature_almudena.iloc[0]
    signature_almudena = signature_almudena.drop("")
    signature_almudena.index.name = "Genes"
    signature_almudena.columns.name = None
    return signature_almudena


def map_hgnc_to_one_ensg(gene_names, adata):
    """
    If a HGNC symbol map to multiple ENSG symbols, choose the one that is in the
    single cell dataset.
    If the HGNC symbol maps to multiple ENSG symbols even inside the scRNAseq dataset,
    then the last one is chosen (no rationale).
    """
    chosen_gene = None
    for gene_name in gene_names:
        if gene_name in adata.var_names:
            chosen_gene = gene_name
    return chosen_gene


def map_hgnc_to_ensg(genes, adata):
    """
    Map the HGNC symbols from the signature matrix to the corresponding ENSG symbols
    of the scRNAseq dataset.
    """
    ensg_names = []
    for gene in genes.index:
        if len(genes.loc[gene].shape) > 1:
            # then one hgnc has multiple ensg lines in the dataframe
            gene_names = genes.loc[gene, "ensembl.gene"]
            gene_name = map_hgnc_to_one_ensg(gene_names, adata)
            if gene_name not in ensg_names:  # for duplicates
                ensg_names.append(gene_name)
        elif genes.loc[gene, "notfound"] is True:
            # then the hgnc gene cannot be mapped to ensg
            ensg_names.append("notfound")
        elif genes.loc[gene, "ensembl.gene"] != genes.loc[gene, "ensembl.gene"]:
            # then one hgnc gene has multiple ensg mappings in one line of the dataframe
            ensembl = genes.loc[gene, "ensembl"]
            gene_names = [ensembl[i]["gene"] for i in range(len(ensembl))]
            gene_name = map_hgnc_to_one_ensg(gene_names, adata)
            ensg_names.append(gene_name)
        else:
            # then one hgnc corresponds to one ensg
            ensg_names.append(genes.loc[gene, "ensembl.gene"])
    return ensg_names


def perform_nnls(signature, averaged_data):
    """Perform deconvolution using the nnls method.
    It will be performed as many times as there are samples in averaged_data.
    """
    deconv = LinearRegression(positive=True).fit(signature, averaged_data.T)
    deconv_results = pd.DataFrame(
        deconv.coef_, index=averaged_data.index, columns=signature.columns
    )
    deconv_results = deconv_results.div(
        deconv_results.sum(axis=1), axis=0
    )  # to sum up to 1
    return deconv_results


def compute_correlations(deconv_results, ground_truth_fractions):
    """Compute n_sample pairwise correlations between the deconvolution results and the
    ground truth fractions.
    """
    deconv_results = deconv_results[ground_truth_fractions.columns] # to align order of columsn
    correlations = [scipy.stats.pearsonr(ground_truth_fractions.iloc[i], 
                                         deconv_results.iloc[i]).statistic 
                                         for i in range(len(deconv_results))]
    correlations = pd.DataFrame({"correlations": correlations})
    correlations["Method"] = "nnls"  # add all the deconv methods
    return correlations
