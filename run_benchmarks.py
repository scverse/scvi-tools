#!/usr/bin/env python

"""Run all the benchmarks with specific parameters"""
import argparse

from scvi.benchmark import run_benchmarks
from scvi.dataset import load_datasets

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--epochs", type=int, default=250, help="how many times to process the dataset")
    parser.add_argument("--dataset", type=str, default="cortex", help="which dataset to process")
    parser.add_argument("--nobatches", action='store_true', help="whether to ignore batches")
    parser.add_argument("--nocuda", action='store_true',
                        help="whether to use cuda (will apply only if cuda is available")
    parser.add_argument("--benchmark", action='store_true',
                        help="whether to use cuda (will apply only if cuda is available")
    parser.add_argument("--url", help="the url for downloading dataset")
    parser.add_argument("--gene_row", help="the row name for genes")
    parser.add_argument("--cluster_col", help="the column name for cell types")
    parser.add_argument("--batch_col", help="the column name for batches")
    args = parser.parse_args()

    dataset = load_datasets(args.dataset, gene_row=args.gene_row,
                            cluster_col=args.cluster_col, batch_col=args.batch_col, url=args.url)
    run_benchmarks(dataset, n_epochs=args.epochs, use_batches=(not args.nobatches), use_cuda=(not args.nocuda),
                   show_batch_mixing=True, benchmark=args.benchmark)
