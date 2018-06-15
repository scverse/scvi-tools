import loompy
from scvi.dataset import load_datasets
#
# ds = loompy.connect("../scVI-data/retina.loom")
# cell_batches = ds.ca['Batch_id']
# labels = ds.ca['Clusters']
#
# genes = ds.ra['Gene']
#
# data = np.array(ds[:, :])  # change matrix to cells by genes
# ds.close()
#
# filename = "retine.loom"
# row_attrs = {'Gene': genes}
# col_attrs = {'ClusterID': labels, 'BatchID': cell_batches}
# loompy.create('../scVI-data/retina.loom', data, row_attrs, col_attrs)

cortex_dataset = load_datasets("Cortex.loom", save_path='data/',
                                   url='http://loom.linnarssonlab.org/clone/Previously%20Published/Cortex.loom')
cortex_dataset.subsample_cells(50)
cortex_dataset.subsample_genes(10)

data = cortex_dataset.X.T  # Genes by Cells
labels = cortex_dataset.labels  # Columns
genes = cortex_dataset.gene_names  # Rows

filename = "Cortex.loom"
save_path = "tests/data/"

row_attrs = {'Gene': genes}
col_attrs = {'ClusterID': labels}
loompy.create(save_path + filename, data, row_attrs, col_attrs)
