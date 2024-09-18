# Counterfactual prediction

Once we have trained a model to predict a variable of interest or a generative model to learn the data distribution, we are often interested in making predictions for new samples. However, predictions over test samples may not reveal exactly what the model has learned about how the input features relate to the target variable of interest. For example, we may want to answer the question: How would the model predict the expression levels of gene Y in cell Z if gene X is knocked out? Even if we do not have an data point corresponding to this scenario, we can instead perturb the input to see what the model reports.

:::{warning}
We are using the term "counterfactual prediction" here rather loosely. In particular, we are not following the rigorous definition of counterfactual prediction in the causality literature[^ref1]. While closely related in spirit, we are making counterfactual queries with statistical models to gain some insight into what the model has learned about the data distribution.
:::

:::{figure} figures/counterfactual_cartoon.svg
:align: center
:alt: Cartoon of the counterfactual prediction task across two conditions.
:class: img-fluid

Cartoon of the counterfactual prediction task across two conditions. This counterfactual prediction can be thought of as an interpolation of nearby points in the feature space originating from condition B.
:::

## Preliminaries

Suppose we have a trained model $f_\theta$ that takes in a data point $x$ (e.g., gene expression counts) and a condition $c$ (e.g., treatment group) and returns a prediction $\hat{y}$.
Each data point takes the form of a tuple $(x,c) \in \mathcal{D}$.
We can define a *counterfactual query* as a pair $(x,c')$ where $c' \neq c$,
and the respective model output as the *counterfactual prediction*, $\hat{y}' = f_\theta(x,c')$.

We separate $c$ here out from $x$ to make the counterfactual portion of the query explicit, but it can be thought of as another dimension of $x$.

## In-distribution vs. out-of-distribution

Since we are working with statistical models rather than causal models, we have to be careful when we can rely on counterfactual predictions. At a high level, if we assume the true function relating the features to the target is smooth, we can trust counterfactual predictions for queries that are similar to points in the training data.

Say we have a counterfactual query $(x,c')$, and we have data points in the training set $(x',c')$ (i.e., $\|x - x'\|$ is small).
If our model predicts the $y$ for $(x', c')$ well,
we can reasonably trust the counterfactual prediction for $(x,c')$.
Otherwise, if $(x,c')$ is very different from any point in the training data
with condition $c'$, we cannot make any guarantees about the accuracy of the counterfactual prediction.
Dimensionality reduction techniques or harmonization methods may help create more overlap between the features $x$ across the conditions, setting the stage for more reliable counterfactual predictions.

## Applications

The most direct application of counterfactual prediction in scvi-tools can be found in the `transform_batch` kwarg of the {func}`~scvi.model.SCVI.get_normalized_expression` function. In this case, we can pass in a counterfactual batch label to get a prediction of what the normalized expression would be for a cell if it were a member of that batch. This can be useful if one wants to compare cells across different batches in the gene space.

The described approach to counterfactual prediction has also been used in a variety of applications, including:
- characterizing cell-type-specific sample-level effects [^ref2]
- predicting chemical perturbation responses in different cell types [^ref2][^ref3]
- predicting infection/perturbation responses across species [^ref4]

For more details on how counterfactual prediction is used in another method implemented in scvi-tools, see the {doc}`/user_guide/models/mrvi`.

[^ref1]:
    Judea Pearl. Causality. Cambridge university press, 2009.
[^ref2]:
     Pierre Boyeau, Justin Hong, Adam Gayoso, Martin Kim, Jose L McFaline-Figueroa, Michael Jordan, Elham Azizi, Can Ergen, Nir Yosef (2024),
    _Deep generative modeling of sample-level heterogeneity in single-cell genomics_,
    [bioRxiv](https://doi.org/10.1101/2022.10.04.510898).
[^ref3]:
    Mohammad Lotfollahi, Anna Klimovskaia Susmelj, Carlo De Donno, Leon Hetzel, Yuge Ji, Ignacio L Ibarra, Sanjay R Srivatsan, Mohsen Naghipourfar, Riza M Daza, Beth Martin, Jay Shendure, Jose L McFaline‐Figueroa, Pierre Boyeau, F Alexander Wolf, Nafissa Yakubova, Stephan Günnemann, Cole Trapnell, David Lopez‐Paz, Fabian J Theis (2023),
    _Predicting cellular responses to complex perturbations in high‐throughput screens_,
    [Molecular Systems Biology](https://doi.org/10.15252/msb.202211517).
[^ref4]:
    Mohammad Lotfollahi, F Alexander Wolf, Fabian J Theis (2019),
    _scGen predicts single-cell perturbation responses_,
    [Nature Methods](https://doi.org/10.1038/s41592-019-0494-8).
