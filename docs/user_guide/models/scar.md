# scAR

**scAR** [^ref1] (single-cell Ambient Remover) is a deep learning model for removal of the ambient signals in droplet-based single cell omics.

This model was ported from another [Github](https://github.com/Novartis/scar).

## Ambient RNA removal

```
>>> adata = anndata.read_h5ad(path_to_anndata)
>>> raw_adata = anndata.read_h5ad(path_to_raw_anndata)
>>> scvi_external.SCAR.setup_anndata(adata, batch_key="batch")
>>> scvi_external.SCAR.get_ambient_profile(adata=adata, raw_adata=raw_adata, prob=0.995)
>>> vae = scvi_external.SCAR(adata, ambient_profile="ambient_profile")
>>> vae.train()
>>> adata.obsm["X_scAR"] = vae.get_latent_representation()
>>> adata.layers['denoised'] = vae.get_denoised_counts()
```

## Estimating the ambient profile

### Option 1: Calculate the ambient profile using the `get_ambient_profile` method

`get_ambient_profile` is inspired by **EmptyDrops** [^ref2], which assumes that cell-free droplets are sampled from a multinomial distribution.

1. It first calculates an initial ambient profile by averaging all droplets in raw_adata
2. It tests whether droplets fit the multinomial distribution (the ambient profile as the prob parameter). The relatively high probabilities suggest cell-free droplets
3. It re-calculates the ambient profile using the identified cell-free droplets
4. It repeats step 2 and step 3 for iterations
5. The final ambient profile is saved in adata.varm

```
>>> adata = anndata.read_h5ad(path_to_anndata)
>>> raw_adata = anndata.read_h5ad(path_to_raw_anndata)
>>> scvi_external.SCAR.get_ambient_profile(adata=adata, raw_adata=raw_adata, prob=0.995)
```

### Option 2: Calculate the ambient profile using a kneeplot

This option is based on total counts of droplets.

1. It first plots a kneeplot
2. It identifies subpopulations of droplets with arbitrary thresholds
3. It calculates the ambient profile from the subpopulation of 'cell-free droplets'

```
>>> all_droplets = pd.DataFrame(raw_adata.X.sum(axis=1), index=raw_adata.obs_names, columns=['total_counts'])
>>> all_droplets['droplets'] = 'cell-free droplets'
>>> all_droplets['droplets'] = all_droplets['droplets'].mask(all_droplets['total_counts']>min_counts, 'other droplets')
>>> all_droplets['droplets'] = all_droplets['droplets'].mask(all_droplets.index.isin(adata.obs_names), 'cells')
>>> all_droplets.index.name = 'barcode'
>>> all_droplets = all_droplets.sort_values(by='total_counts', ascending=False).reset_index().rename_axis("rank").reset_index()
>>> all_droplets = all_droplets.loc[all_droplets['total_counts']>0]
>>> all_droplets = all_droplets.set_index('barcode').rename_axis('cells')

>>> plt.figure(figsize=(3, 2), dpi=150)
>>> ax = sns.lineplot(
>>>     data=all_droplets,
>>>     x='rank',
>>>     y='total_counts',
>>>     hue='droplets',
>>>     hue_order=['other droplets', 'cell-free droplets', 'cells'],
>>>     palette=sns.color_palette()[-3:],
>>>     markers=False,
>>>     lw=2,
>>>     ci=None
>>> )
>>> ax.set_xscale('log')
>>> ax.set_yscale('log')
>>> ax.set_xlabel('sorted droplets');
>>> ax.legend(loc='lower left', ncol=1, title=None, frameon=False)
>>> ax.set_title(f'kneeplot')
>>> sns.set_palette('muted')
>>> sns.set_style('ticks')
>>> sns.despine(offset=10, trim=False)

>>> cell_free = raw_adata[raw_adata.obs_names.isin(all_droplets[all_droplets['droplets']=='cell-free droplets'].index)].copy()
>>> cell_free = cell_free[:, adata.var_names]
>>> # average and normalize the transcript in cell-free droplets
>>> ambient_profile = pd.DataFrame((cell_free.X.sum(axis=0)/cell_free.X.sum()).A1, index=adata.var_names, columns=['ambient profile'])
```

[^ref1]:
    Caibin Sheng, Rui Lopes, Gang Li, Sven Schuierer, Annick Waldt, Rachel Cuttat, Slavica Dimitrieva, Audrey Kauffmann, Eric Durand, Giorgio G. Galli, Guglielmo Roma, Antoine de Weck (2022),
    _Probabilistic machine learning ensures accurate ambient denoising in droplet-based single-cell omics_,
    [bioRxiv](https://www.biorxiv.org/content/10.1101/2022.01.14.476312v4).

[^ref2]:
    Aaron T. L. Lun, Samantha Riesenfeld, Tallulah Andrews, The Phuong Dao, Tomas Gomes, participants in the 1st Human Cell Atlas Jamboree & John C. Marioni (2019),
    _EmptyDrops: distinguishing cells from empty droplets in droplet-based single-cell RNA sequencing data_,
    [biomedcentral](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1662-y)
