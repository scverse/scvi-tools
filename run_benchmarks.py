#!/usr/bin/env python

"""Run all the benchmarks with specific parameters"""
import argparse

from scvi.benchmark import run_benchmarks
from scvi.dataset import BrainLargeDataset, CortexDataset, SyntheticDataset, \
    RetinaDataset, CbmcDataset, BrainSmallDataset, HematoDataset, PbmcDataset, LoomDataset


def load_datasets(dataset_name, save_path='data/', url=None):
    if dataset_name == 'synthetic':
        gene_dataset = SyntheticDataset()
    elif dataset_name == 'cortex':
        gene_dataset = CortexDataset()
    elif dataset_name == 'brain_large':
        gene_dataset = BrainLargeDataset(subsample_size=128, save_path=save_path)
    elif dataset_name == 'retina':
        gene_dataset = RetinaDataset(save_path=save_path)
    elif dataset_name == 'cbmc':
        gene_dataset = CbmcDataset(save_path=save_path)
    elif dataset_name == 'brain_small':
        gene_dataset = BrainSmallDataset(save_path=save_path)
    elif dataset_name == 'hemato':
        gene_dataset = HematoDataset(save_path=save_path)
    elif dataset_name == 'pbmc':
        gene_dataset = PbmcDataset(save_path=save_path)
    elif dataset_name[-5:] == ".loom":
        gene_dataset = LoomDataset(filename=dataset_name, save_path=save_path, url=url)
    else:
        raise "No such dataset available"
    return gene_dataset


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--epochs", type=int, default=250, help="how many times to process the dataset")
    parser.add_argument("--dataset", type=str, default="cortex", help="which dataset to process")
    parser.add_argument("--nobatches", action='store_true', help="whether to ignore batches")
    parser.add_argument("--nocuda", action='store_true',
                        help="whether to use cuda (will apply only if cuda is available")
    parser.add_argument("--benchmark", action='store_true',
                        help="whether to use cuda (will apply only if cuda is available")
    parser.add_argument("--url", type=str, help="the url for downloading dataset")
    args = parser.parse_args()

    dataset = load_datasets(args.dataset, url=args.url)
    run_benchmarks(dataset, n_epochs=args.epochs, use_batches=(not args.nobatches), use_cuda=(not args.nocuda),
                   show_batch_mixing=True, benchmark=args.benchmark)
