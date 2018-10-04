import pandas as pd
import numpy as np
res = pd.read_table('../PopRemove/res.txt')
vae = res[res['model_type']=='vae']
seurat = res[res['model_type']=='Seurat']
