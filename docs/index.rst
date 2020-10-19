========================
scvi-tools documentation
========================

scvi-tools (single-cell variational inference tools) is a package for end-to-end analysis of single-cell omics data. The package is primarily developed and maintained by the `Yosef Lab
<https://yoseflab.github.io/>`_ at UC Berkeley and is composed of several deep generative models for omics data analysis:

* scVI for analysis of single-cell RNA-seq data [Lopez18]_
* scANVI for cell annotation of scRNA-seq data using semi-labeled examples [Xu19]_
* totalVI for analysis of CITE-seq data [GayosoSteier20]_
* gimVI for imputation of missing genes in spatial transcriptomics from scRNA-seq data [Lopez19]_
* AutoZI for assessing gene-specific levels of zero-inflation in scRNA-seq data [Clivio19]_
* LDVAE for an interpretable linear factor model version of scVI [Svensson20]_

These models are able to simultaneously perform many downstream tasks such as learning low-dimensional cell representations, harmonizing datasets from different experiments, and identifying differential expressed features [Boyeau19]_. By levaraging advances in stochastic optimization, these models scale to millions of cells.

* If you find a model useful for your research, please consider citing the corresponding publication.

.. important:: ``scvi`` is now ``scvi-tools``.
   If you'd like to view documentation for ``scvi``, please change the documentation version using the menu at the bottom right (versions <= 0.6.8).

.. raw:: html

    <div class="container">
        <div class="row">
            <div class="col-lg-6 col-md-6 col-sm-6 col-xs-12 d-flex">
                <div class="card text-center intro-card shadow">
                <img src="_static/computer-24px.svg" class="card-img-top" alt="getting started with scvi action icon" height="52">
                <div class="card-body flex-fill">
                    <h5 class="card-title">Getting started</h5>
                    <p class="card-text">New to <em>scvi-tools</em>? Check out the installation guide.
                    </p>

.. container:: custom-button

    :doc:`To the installation guide<installation>`

.. raw:: html

                </div>
                </div>
            </div>
            <div class="col-lg-6 col-md-6 col-sm-6 col-xs-12 d-flex">
                <div class="card text-center intro-card shadow">
                <img src="_static/play_circle_outline-24px.svg" class="card-img-top" alt="scvi user guide action icon" height="52">
                <div class="card-body flex-fill">
                    <h5 class="card-title">User guide</h5>
                    <p class="card-text">The tutorials provide in-depth information on running scvi-tools models.</p>

.. container:: custom-button

    :doc:`To the user guide<user_guide/index>`

.. raw:: html

                </div>
                </div>
            </div>
            <div class="col-lg-6 col-md-6 col-sm-6 col-xs-12 d-flex">
                <div class="card text-center intro-card shadow">
                <img src="_static/library_books-24px.svg" class="card-img-top" alt="api of scvi action icon" height="52">
                <div class="card-body flex-fill">
                    <h5 class="card-title">API reference</h5>
                    <p class="card-text">The API reference contains a detailed description of
                    the scvi-tools API.</p>

.. container:: custom-button

    :doc:`To the API reference<api/index>`

.. raw:: html

                </div>
                </div>
            </div>
            <div class="col-lg-6 col-md-6 col-sm-6 col-xs-12 d-flex">
                <div class="card text-center intro-card shadow">
                <img src="_static/code-24px.svg" class="card-img-top" alt="contribute to pandas action icon" height="52">
                <div class="card-body flex-fill">
                    <h5 class="card-title">Developer guide</h5>
                    <p class="card-text">Want to improve scvi-tools? The contributing guidelines
                    will guide you through the process.</p>

.. container:: custom-button

    :doc:`To the development guide<contributing>`

.. raw:: html

                </div>
                </div>
            </div>
        </div>
    </div>


.. toctree::
   :maxdepth: 3
   :hidden:
   :titlesonly:

   installation
   user_guide/index
   api/index
   contributing
   history
   authors
   references
   Discussion <https://discourse.scvi-tools.org>
