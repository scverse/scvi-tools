#!/usr/bin/env python


"""Run all the benchmarks with specific parameters"""
import argparse
import time

from scvi.benchmark import run_benchmarks
from scvi.dataset import load_datasets

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--epochs", type=int, default=250, help="how many times to process the dataset")
    parser.add_argument("--dataset", type=str, default="cortex", help="which dataset to process")
    parser.add_argument("--nobatches", action='store_true', help="whether to ignore batches")

    args = parser.parse_args()
    gene_dataset_train, gene_dataset_test = load_datasets(args.dataset)
    start = time.time()
    run_benchmarks(gene_dataset_train, gene_dataset_test, n_epochs=args.epochs,
                   use_batches=(not args.nobatches))
    end = time.time()
    print("Total runtime for " + str(args.epochs) + " epochs is: " + str((end - start))
          + " seconds for a mean per epoch runtime of " + str((end - start) / args.epochs) + " seconds.")
