import warnings
from collections import namedtuple

import numpy as np
import rpy2.robjects as ro
import rpy2.robjects.numpy2ri
from rpy2.rinterface import RRuntimeWarning

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


    def predict_scmap_cell(self, data_test, labels_test):
        if hasattr(data_test, 'A'):
            data_train = data_test.A
        self.projection = 'test'
        self.create_sce_object(data_test, np.arange(data_test.shape[1]).astype(np.str),
                               labels_test, self.projection)
        ro.r("sce <- indexCell(%s)" % self.reference)
        ro.r("result<-scmapCell(%s, list(metadata(%s)$scmap_cell_index))"
             % (self.projection, self.reference))

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

