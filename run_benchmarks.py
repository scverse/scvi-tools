#!/usr/bin/env python

"""Run all the benchmarks with specific parameters"""
import argparse

from scvi.benchmark import harmonization_benchmarks, \
    annotation_benchmarks, all_benchmarks
from scvi.dataset import BrainLargeDataset, CortexDataset, SyntheticDataset, CsvDataset, \
    RetinaDataset, BrainSmallDataset, HematoDataset, LoomDataset, AnnDataset, CbmcDataset, PbmcDataset
from scvi.inference import UnsupervisedTrainer, SemiSupervisedTrainer
from scvi.models import VAE, VAEC, SCANVI


def load_datasets(dataset_name, save_path='data/', url=None):
    if dataset_name == 'synthetic':
        gene_dataset = SyntheticDataset()
    elif dataset_name == 'cortex':
        gene_dataset = CortexDataset()
    elif dataset_name == 'brain_large':
        gene_dataset = BrainLargeDataset(save_path=save_path)
    elif dataset_name == 'retina':
        gene_dataset = RetinaDataset(save_path=save_path)
    elif dataset_name == 'cbmc':
        gene_dataset = CbmcDataset(save_path=save_path)
    elif dataset_name == 'brain_small':
        gene_dataset = BrainSmallDataset(save_path=save_path)
    elif dataset_name == 'hemato':
        gene_dataset = HematoDataset(save_path='data/HEMATO/')
    elif dataset_name == 'pbmc':
        gene_dataset = PbmcDataset(save_path=save_path)
    elif dataset_name[-5:] == ".loom":
        gene_dataset = LoomDataset(filename=dataset_name, save_path=save_path, url=url)
    elif dataset_name[-5:] == ".h5ad":
        gene_dataset = AnnDataset(dataset_name, save_path=save_path, url=url)
    elif ".csv" in dataset_name:
        gene_dataset = CsvDataset(dataset_name, save_path=save_path)
    else:
        raise "No such dataset available"
    return gene_dataset


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


available_models = {
    'VAE': VAE,
    'VAEC': VAEC,
    'SVAEC': SCANVI
}

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--epochs", type=int, default=250, help="how many times to process the gene_dataset")
    parser.add_argument("--dataset", type=str, default="cortex", help="which gene_dataset to process")
    parser.add_argument("--model", type=str, default="VAE", help="the model to use")
    parser.add_argument("--nobatches", action='store_true', help="whether to ignore batches")
    parser.add_argument("--nocuda", action='store_true',
                        help="whether to use cuda (will apply only if cuda is available")
    parser.add_argument("--harmonization", action='store_true',
                        help="whether to use cuda (will apply only if cuda is available")
    parser.add_argument("--annotation", action='store_true',
                        help="whether to use cuda (will apply only if cuda is available")
    parser.add_argument("--all", action='store_true',
                        help="whether to use cuda (will apply only if cuda is available")
    parser.add_argument("--benchmark", action='store_true',
                        help="whether to use cuda (will apply only if cuda is available")
    parser.add_argument("--url", type=str, help="the url for downloading gene_dataset")
    args = parser.parse_args()

    n_epochs = args.epochs
    use_cuda = not args.nocuda
    if args.all:
        all_benchmarks(n_epochs=n_epochs, use_cuda=use_cuda)
    elif args.harmonization:
        harmonization_benchmarks(n_epochs=n_epochs, use_cuda=use_cuda)
    elif args.annotation:
        annotation_benchmarks(n_epochs=n_epochs, use_cuda=use_cuda)
    else:
        dataset = load_datasets(args.dataset, url=args.url)
        model = available_models[args.model](dataset.nb_genes, dataset.n_batches*args.nobatches, dataset.n_labels)
        trainer_cls = UnsupervisedTrainer if args.model == 'VAE' else SemiSupervisedTrainer
        trainer = trainer_cls(model, dataset, use_cuda=use_cuda)
        trainer.train(n_epochs=n_epochs)
