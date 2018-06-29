from . import BrainLargeDataset


def benchmark_hyperparameters(gene_dataset):
    """
    Possibly specify hyper parameters according to Table2 in "Bayesian Inference for a Generative Model of
    transcriptome Profiles from Single-cell RNA Sequencing"
    :param gene_dataset:
    :return:
    """
    hyperparameters = dict()
    if isinstance(gene_dataset, BrainLargeDataset):
        if gene_dataset.nb_genes > 10000 and gene_dataset.nb_genes < 15000:
            hyperparameters["n_layers"] = 2
        elif gene_dataset.nb_genes > 15000:
            hyperparameters["n_layers"] = 3
            hyperparameters["n_hidden"] = 256

    return hyperparameters
