import os
import urllib.request

import numpy as np
import pandas as pd

from .dataset import GeneExpressionDataset
from .utils import train_test_split


class CortexDataset(GeneExpressionDataset):
    def __init__(self):
        # Generating samples according to a ZINB process
        self.save_path = 'data/'
        self.download_name = 'expression.bin'
        self.final_name = 'expression_train.npy'
        self.download_and_preprocess()
        super(CortexDataset, self).__init__([np.load(self.save_path + self.final_name)])

    def download(self):
        url = "https://storage.googleapis.com/linnarsson-lab-www-blobs/blobs/cortex/expression_mRNA_17-Aug-2014.txt"
        r = urllib.request.urlopen(url)
        # total_size = int(r.headers['content-length']) / 1000
        print("Downloading Cortex data")

        def readIter(f, blocksize=1000):
            """Given a file 'f', returns an iterator that returns bytes of
            size 'blocksize' from the file, using read()."""
            while True:
                data = f.read(blocksize)
                if not data:
                    break
                yield data

        # Create the path to save the data

        if not os.path.exists(self.save_path):
            os.makedirs(self.save_path)

        with open(self.save_path + self.download_name, 'wb') as f:
            for data in readIter(r):  # tqdm(readIter(r), total=total_size, unit='KB', unit_scale=False):
                f.write(data)

    def preprocess(self):
        print("Preprocessing Cortex data")
        X = pd.read_csv(self.save_path + self.download_name, sep="\t", low_memory=False).T
        clusters = np.array(X[7], dtype=str)[2:]
        cell_types, labels = np.unique(clusters, return_inverse=True)
        gene_names = np.array(X.iloc[0], dtype=str)[10:]
        X = X.loc[:, 10:]
        X = X.drop(X.index[0])
        expression_data = np.array(X, dtype=np.int)[1:]

        # keep the most variable genes according to the Biscuit ICML paper
        selected = np.std(expression_data, axis=0).argsort()[-558:][::-1]
        expression_data = expression_data[:, selected]
        gene_names = gene_names[selected].astype(str)

        # train test split for log-likelihood scores
        expression_train, expression_test, c_train, c_test = train_test_split(expression_data, labels)

        np.save(self.save_path + self.final_name, expression_train)
        np.save(self.save_path + 'expression_test.npy', expression_test)

    def download_and_preprocess(self):
        if not os.path.exists(self.save_path + self.final_name):
            if not os.path.exists(self.save_path + self.download_name):
                self.download()
            self.preprocess()
