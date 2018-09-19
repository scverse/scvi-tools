import warnings
from collections import namedtuple

import numpy as np
import rpy2.robjects as ro
import rpy2.robjects.numpy2ri
from rpy2.rinterface import RRuntimeWarning

from scvi.dataset import GeneExpressionDataset, SemiSupervisedDataLoaders
from scvi.metrics.classification import compute_accuracy_tuple

Accuracy = namedtuple('Accuracy',
                      ['accuracy',
                       'unclassified_rate',
                       'accuracy_over_known_classes',
                       'accuracy_over_known_classes_and_assigned'])


def convert_labels_str(str_labels):
    str_labels[str_labels == 'unassigned'] = -1
    return str_labels.astype(np.int)


def convert_labels_levels(r_indices, levels):
    levels = (levels.astype(np.float32)).astype(np.int)
    return np.array([int(levels[int(l) - 1]) if l != "unassigned" else -1 for l in r_indices.astype(np.int)])


class SCMAP():
    # An algorithm for annotation
    def __init__(self):
        warnings.filterwarnings("ignore", category=RRuntimeWarning)
        rpy2.robjects.numpy2ri.activate()
        ro.r["library"]("scmap")
        ro.r["library"]("SingleCellExperiment")
        ro.r["library"]("matrixStats")

        self.n_features=100
        self.threshold =0

    def set_parameters(self,n_features=100, threshold=0):
        self.n_features = n_features
        self.threshold = threshold

    def create_dataset(self, path):
        print("Reading rds")
        ro.r("sce<-readRDS('%s')" % path)
        print("Extracting log counts")
        log_counts = ro.r("logcounts(sce)")
        print("Transforming log count to counts")
        counts = (np.exp(log_counts * np.log(2)) - 1).T.astype(np.int)
        gene_symbols = ro.r("rowData(sce)$feature_symbol")
        labels = ro.r("colData(sce)$cell_type1")
        labels_levels = ro.r("levels(colData(sce)$cell_type1)")
        if labels_levels is not rpy2.rinterface.NULL:
            labels = np.array([labels_levels[int(l) - 1] for l in labels])

        cell_types = list(np.unique(labels))
        labels = np.array([cell_types.index(l) for l in labels])

        valid_idx = (counts.sum(axis=1) > 10).ravel()  # Filter bad quality cells
        counts = counts[valid_idx]
        labels = labels[valid_idx]
        gene_expression_dataset = GeneExpressionDataset(
            *GeneExpressionDataset.get_attributes_from_matrix(counts, labels=labels), cell_types=cell_types
        )
        gene_expression_dataset.gene_symbols = gene_symbols
        return gene_expression_dataset

    def create_sce_object(self, matrix, gene_names, labels, name):
        n_samples, nb_genes = matrix.shape
        r_matrix = ro.r.matrix(matrix.T, nrow=nb_genes, ncol=n_samples)
        ro.r.assign("counts", r_matrix)
        ro.r.assign("gene_names", ro.StrVector(gene_names))
        ro.r("counts<-as.data.frame(counts, row.names=gene_names)")

        ro.r.assign("barcodes_cells", ro.StrVector(["cell_" + str(i) for i in range(n_samples)]))
        ro.r("colnames(counts)<-barcodes_cells")

        if labels is not None:
            ro.r.assign("labels", ro.StrVector(labels))
            ro.r("barcodes_cells<-as.data.frame(labels, row.names=barcodes_cells, col.names=c('cell_type1'))")
            ro.r("colnames(barcodes_cells)<-c('cell_type1')")
            ro.r("%s <- SingleCellExperiment(assays=list(counts=as.matrix(counts)), colData=barcodes_cells)" % name)
        else:
            ro.r("%s <- SingleCellExperiment(assays=list(counts=as.matrix(counts)))" % name)
        ro.r("rowData(%s)$feature_symbol<-rownames(%s)" % (name, name))  # For any new custom dataset.
        ro.r("logcounts(%s) <- log2(counts(%s) + 1)" % (name, name))
        assert np.sum(self.get_labels(name) != labels) == 0, "Labels not matching"  # Sanity check R/python
        print("SCE object : %s created" % name)

    def select_features(self, reference, n_features=500):
        ro.r("%s<-selectFeatures(%s,  n_features=%d)" % (reference, reference, n_features))
        scmap_features = ro.r("rowData(%s)$scmap_features" % reference)
        print("%i/%i features selected" % (np.sum(scmap_features), len(scmap_features)))

    def scmap_cluster(self, reference, projection, threshold=0, n_features=500):
        self.select_features(reference, n_features=n_features)

        ro.r("%s<-indexCluster(%s)" % (reference, reference))
        ro.r("result<-scmapCluster(%s, list(metadata(%s)$scmap_cluster_index), threshold=%.1f)"
             % (projection, reference, threshold))  # list(metadata(sce_reference)$scmap_cluster_index))")

        self.probs = ro.r("result$scmap_cluster_siml")

        self.labels_pred = convert_labels_str(ro.r("result$scmap_cluster_labs"))  # 'unassigned' are included
        self.combined_labels_pred = convert_labels_str(ro.r("result$combined_labs"))

        return self.labels_pred

    def fit_scmap_cluster(self, data_train, labels_train):
        if hasattr(data_train, 'A'):
            data_train = data_train.A
        self.reference = 'train'
        self.create_sce_object(data_train, np.arange(data_train.shape[1]).astype(np.str),
                               labels_train.astype(np.int), self.reference)

        self.select_features(self.reference, n_features=self.n_features)

        ro.r("%s<-indexCluster(%s)" % (self.reference, self.reference))

    def predict_scmap_cluster(self, data_test, labels_test):
        if hasattr(data_test, 'A'):
            data_train = data_test.A
        self.projection = 'test'
        self.create_sce_object(data_test, np.arange(data_test.shape[1]).astype(np.str),
                               labels_test, self.projection)
        ro.r("result<-scmapCluster(%s, list(metadata(%s)$scmap_cluster_index), threshold=%.1f)"
             % (self.projection, self.reference, self.threshold))

        self.probs = ro.r("result$scmap_cluster_siml")

        self.labels_pred = convert_labels_str(ro.r("result$scmap_cluster_labs"))  # 'unassigned' are included
        self.combined_labels_pred = convert_labels_str(ro.r("result$combined_labs"))
        return self.labels_pred

    def score(self, data_test, labels_test):
        labels_pred = self.predict_scmap_cluster(data_test, labels_test)
        self.test_tuple = compute_accuracy_tuple(labels_test, labels_pred)
        return np.mean(labels_pred == labels_test)

    def scmap_cell(self, reference, projection, w=10, threshold=0, n_features=500):
        self.select_features(reference, n_features=n_features)

        ro.r("%s<-indexCell(%s)" % (reference, reference))
        ro.r("scmapCell_result<-scmapCell(%s, list(metadata(%s)$scmap_cell_index),w=%d)" % (projection, reference, w))

        ro.r("result<-scmapCell2Cluster(scmapCell_result, cluster_list=list(colData(%s)$cell_type1, "
             "w=%d, threshold=%.1f))" % (reference, w, threshold))

        self.probs = ro.r("result$scmap_cluster_siml")

        self.levels_reference = ro.r('levels(colData(%s)$cell_type1)' % reference)

        self.labels_pred = convert_labels_levels(ro.r("result$scmap_cluster_labs"), self.levels_reference)
        self.combined_labels_pred = convert_labels_levels(ro.r("result$combined_labs"), self.levels_reference)
        return self.labels_pred

    def scmap_cluster_table(self, d1, d2, d12, n_labelled_samples_per_class=10, threshold=0, n_features=500):
        r"""
        Given two different gene expression dataset
        """
        data_loaders_1 = SemiSupervisedDataLoaders(d1, n_labelled_samples_per_class=n_labelled_samples_per_class)
        (X_train_1, labels_train_1), = data_loaders_1.raw_data(data_loaders_1['labelled'])
        (X_test_1, labels_test_1), = data_loaders_1.raw_data(data_loaders_1['sequential'])
        self.create_sce_object(X_train_1, d1.gene_names, labels_train_1, 'd1_train')
        self.create_sce_object(X_test_1, d1.gene_names, labels_test_1, 'd1_test')

        data_loaders_2 = SemiSupervisedDataLoaders(d2,
                                                   n_labelled_samples_per_class=n_labelled_samples_per_class)
        (X_train_2, labels_train_2), = data_loaders_2.raw_data(data_loaders_2['labelled'])
        (X_test_2, labels_test_2), = data_loaders_2.raw_data(data_loaders_2['sequential'])
        self.create_sce_object(X_train_2, d2.gene_names, labels_train_2, 'd2_train')
        self.create_sce_object(X_test_2, d2.gene_names, labels_test_2, 'd2_test')

        data_loaders_12 = SemiSupervisedDataLoaders(d12,
                                                    n_labelled_samples_per_class=n_labelled_samples_per_class)
        (X_train_12, labels_train_12), = data_loaders_12.raw_data(data_loaders_12['labelled'])
        (X_test_12, labels_test_12), = data_loaders_12.raw_data(data_loaders_12['sequential'])
        self.create_sce_object(X_train_12, d2.gene_names, labels_train_12, 'd12_train')
        self.create_sce_object(X_test_12, d2.gene_names, labels_test_12, 'd12_test')
        accuracy_results = np.zeros((3, 3))

        for i, reference in enumerate(['d1', 'd2', 'd12']):
            for j, target in enumerate(['d1', 'd2', 'd12']):
                labels_pred = self.scmap_cluster(
                    reference + '_train', target + '_test', threshold=threshold, n_features=n_features
                )
                labels = self.get_labels(target + '_test')
                accuracy_results[i, j] = np.mean(labels_pred == labels)
                print("Accuracy results :", accuracy_results[i, j])
        print(accuracy_results)
        return accuracy_results

    def scmap_cluster_self(self, d, n_labelled_samples_per_class=10, threshold=0, n_features=500):
        r'''
        Self-classification on an scmap cluster
        :param d:
        :return:
        '''
        data_loaders = SemiSupervisedDataLoaders(d, n_labelled_samples_per_class=n_labelled_samples_per_class)
        (X_train, labels_train), = data_loaders.raw_data(data_loaders['labelled'])
        (X_test, labels_test), = data_loaders.raw_data(data_loaders['sequential'])
        self.create_sce_object(X_train, d.gene_names, labels_train, 'd_train')
        self.create_sce_object(X_test, d.gene_names, labels_test, 'd_test')
        labels_pred = self.scmap_cluster(
            'd_train', 'd_test', threshold=threshold, n_features=n_features
        )

        accuracy = np.mean(labels_pred == labels_test)
        print("Accuracy result :", accuracy)
        return accuracy

    def get_levels(self, name):  # TODO : THIS IS VERY FALSE
        return ro.r('levels(colData(%s)$cell_type1)' % name)

    def get_labels(self, name):  # TODO : THIS IS VERY FALSE
        levels = self.get_levels(name)
        r_indices = ro.r("colData(%s)$cell_type1" % name)
        return convert_labels_levels(r_indices, levels)

    def accuracy_tuple(self, y, y_pred, y_train=None):
        if y_train:
            unique_train = np.unique(y_train)
            idx = (sum([(y == i) for i in unique_train])).astype(np.bool)
            idx_strict = np.logical_and(idx, (y_pred != -1))
        a1 = np.nan if not y_train else np.mean(y_pred[idx] == y[idx])
        a2 = np.nan if not y_train else np.mean(y_pred[idx_strict] == y[idx_strict])
        return Accuracy(
            accuracy=np.mean(y_pred == y),
            unclassified_rate=np.mean(y_pred == -1),
            accuracy_over_known_classes=a1,
            accuracy_over_known_classes_and_assigned=a2
        )

    @staticmethod
    def load_rds_file(name, filename):
        ro.r("%s <- readRDS('%s')" % (name, filename))


