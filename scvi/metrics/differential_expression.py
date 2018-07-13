
@no_grad()
@eval_modules()
def most_expressed_genes(vae, data_loader, M_sampling=100, M_permutation=100000, permutation=False):
    r"""
    :param vae: a VAE model object
    :param data_loader: a ``data_loader`` object
    :param M_sampling: number of samples
    :param M_permutation: M_permutation: 10000 - default value in Romain's code
    :param permutation:
    :return: a table of most expressed gene for each cell. Entry[i][j] represents the cell j's gene_expression
    level for gene i. Gene i is the most expressed gene for cell i.
    """
    gene_names = data_loader.dataset.gene_names
    px_scale, all_labels = de_stats(vae, data_loader, M_sampling)

    most_expressed_genes = []
    results = []

    n_cells = vae.n_labels
    for cell_idx in range(n_cells):
        sample_rate_a = (px_scale[all_labels.view(-1) == cell_idx].view(-1, px_scale.size(1))
                         .cpu().detach().numpy())
        sample_rate_b = (px_scale[all_labels.view(-1) != cell_idx].view(-1, px_scale.size(1))
                         .cpu().detach().numpy())
        samples = np.vstack((sample_rate_a, sample_rate_b))

        list_1 = list(np.arange(sample_rate_a.shape[0]))
        list_2 = list(sample_rate_a.shape[0] + np.arange(sample_rate_b.shape[0]))
        if not permutation:
            u, v = np.random.choice(list_1, size=M_permutation), np.random.choice(list_2, size=M_permutation)
        else:
            u, v = (np.random.choice(list_1 + list_2, size=M_permutation),
                    np.random.choice(list_1 + list_2, size=M_permutation))

        first_set = samples[u]
        second_set = samples[v]
        res = np.mean(first_set >= second_set, 0)
        res = np.log(res + 1e-8) - np.log(1 - res + 1e-8)
        results.append(res)
        most_expressed_gene = gene_names[np.argmax(res)]
        most_expressed_genes.append(most_expressed_gene)

    results = [[res[np.where(gene_names == gene_name)[0]][0] for gene_name in most_expressed_genes] for res in results]
    return most_expressed_genes, results
