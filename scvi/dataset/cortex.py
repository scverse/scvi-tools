import csv
import os

import numpy as np
import scipy.sparse as sp_sparse
import urllib.request

from .dataset import GeneExpressionDataset


class CortexDataset(GeneExpressionDataset):
    save_path = 'data/'
    download_name = 'expression.bin'
    data_filename = 'expression_%s.npy'
    labels_filename = 'labels_%s.npy'
    gene_names = 'genes_names.npy'

    def download(self):
        url = "https://storage.googleapis.com/linnarsson-lab-www-blobs/blobs/cortex/expression_mRNA_17-Aug-2014.txt"
        r = urllib.request.urlopen(url)
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
        rows = []
        gene_names = []
        with open(self.save_path + self.download_name, 'r') as csvfile:
            data_reader = csv.reader(csvfile, delimiter='\t')
            clusters = None
            for i, row in enumerate(data_reader):
                if i == 8:  # 7 + 1 in pandas
                    clusters = np.array(row, dtype=str)[2:]
                if i >= 11:  # 10 + 1 in pandas
                    rows.append(row[1:])
                    gene_names.append(row[0])

        cell_types, labels = np.unique(clusters, return_inverse=True)

        expression_data = np.array(rows, dtype=np.int).T[1:]
        gene_names = np.array(gene_names, dtype=np.str)

        selected = np.std(expression_data, axis=0).argsort()[-558:][::-1]
        expression_data = expression_data[:, selected]
        gene_names = gene_names[selected]

        # train test split for log-likelihood scores
        expression_train, expression_test, c_train, c_test = GeneExpressionDataset.train_test_split(expression_data,
                                                                                                    labels)
        np.save(self.save_path + 'expression_train.npy', expression_train)
        np.save(self.save_path + 'expression_test.npy', expression_test)
        np.save(self.save_path + 'labels_train.npy', c_train)
        np.save(self.save_path + 'labels_test.npy', c_test)
        np.save(self.save_path + self.gene_names, gene_names)

    @classmethod
    def download_and_preprocess(self, type='train'):
        if not (os.path.exists(self.save_path + self.data_filename % type) and
                os.path.exists(self.save_path + self.labels_filename % type) and
                os.path.exists(self.save_path + self.gene_names)):
            if not os.path.exists(self.save_path + self.download_name):
                self.download()
            self.preprocess()

    @classmethod
    def get_attributes(self, type='train'):
        return sp_sparse.csr_matrix(np.load(self.save_path + self.data_filename % type)), \
               np.load(self.save_path + self.labels_filename % type), \
               np.load(self.save_path + self.gene_names)

    @classmethod
    def get_dataset(cls, type='train'):
        cls.download_and_preprocess(type=type)
        return cls.from_matrix(*cls.get_attributes(type=type))
