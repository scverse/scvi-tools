User guide
==========

The easiest way to get familiar with scvi-tools is to follow along with our tutorials.
Many are also designed to work seamlessly in Google Colab, a free cloud computing platform. These tutorials have a Colab badge in their introduction. In general, these tutorials are designed to work with the latest installable version of scvi-tools.

To download the tutorials:

1. Follow the Colab link within the tutorial (if available).
2. Download it with the option in the file menu.

.. note:: For questions about using scvi-tools, or broader questions about modeling data, please use our forum_.

.. _forum: https://discourse.scvi-tools.org/

Quick start
-----------

.. raw:: html

    <div class="container">
        <div class="row">
            <div class="col-lg-6 col-md-6 col-sm-6 col-xs-12 d-flex">
                <div class="card text-center intro-card shadow">
                <img src="../_static/tutorials/overview.svg" class="card-img-top" alt="scvi basic tutorial action icon" height="225">
                <div class="card-body flex-fill">
                    <h5 class="card-title">Quick guide to scvi-tools</h5>
                    <p class="card-text">Rapidly learn the basics to run any of the scvi-tools models.
                    </p>

.. container:: custom-button

    :doc:`GO<notebooks/api_overview>`


.. raw:: html

                </div>
                </div>
            </div>
            <div class="col-lg-6 col-md-6 col-sm-6 col-xs-12 d-flex">
                <div class="card text-center intro-card shadow">
                <img src="../_static/tutorials/anndata.svg" class="card-img-top" alt="data loading tutorial action icon" height="225">
                <div class="card-body flex-fill">
                    <h5 class="card-title">Data loading and prep</h5>
                    <p class="card-text">How do I get my data prepared for scvi-tools?</p>

.. container:: custom-button

    :doc:`GO<notebooks/data_loading>`


.. raw:: html

                </div>
                </div>
            </div>
        </div>
    </div>

Tutorials
---------

.. raw:: html

    <div class="container">
    <div id="accordion" class="shadow tutorial-accordion">

        <div class="card tutorial-card">
            <div class="card-header collapsed card-link" data-toggle="collapse" data-target="#collapseOne">
                <div class="d-flex flex-row tutorial-card-header-1">
                    <div class="d-flex flex-row tutorial-card-header-2">
                        <button class="btn btn-dark btn-sm"></button>
                        Analysis of paired RNA and protein data with totalVI
                    </div>
                    <span class="badge gs-badge-link">

:doc:`Straight to tutorial...<notebooks/totalVI>`

.. raw:: html

                    </span>
                </div>
            </div>
            <div id="collapseOne" class="collapse" data-parent="#accordion">
                <div class="card-body">

This is a walkthrough of a totalVI-based analysis pipeline, from dimension reduction to differential expression.

.. image:: ../_static/tutorials/totalvi_cell.svg
   :align: center
   :height: 300px

.. raw:: html

                    <div class="d-flex flex-row">
                        <span class="badge gs-badge-link">

:doc:`To the tutorial <notebooks/totalVI>`

.. raw:: html

                        </span>
                    </div>
                </div>
            </div>
        </div>

        <div class="card tutorial-card">
            <div class="card-header collapsed card-link" data-toggle="collapse" data-target="#collapseTwo">
                <div class="d-flex flex-row tutorial-card-header-1">
                    <div class="d-flex flex-row tutorial-card-header-2">
                        <button class="btn btn-dark btn-sm"></button>
                        Atlas-level integration and label transfer
                    </div>
                    <span class="badge gs-badge-link">

:doc:`Straight to tutorial...<notebooks/harmonization>`

.. raw:: html

                    </span>
                </div>
            </div>
            <div id="collapseTwo" class="collapse" data-parent="#accordion">
                <div class="card-body">

Here we describe how to use scVI and scANVI for integrating data from Tabula Muris.

.. image:: ../_static/tutorials/scanvi.png
   :align: center
   :height: 300px

.. raw:: html

                    <div class="d-flex flex-row">
                        <span class="badge gs-badge-link">

:doc:`To the tutorial <notebooks/harmonization>`


.. raw:: html

                        </span>
                    </div>
                </div>
            </div>
        </div>

        <div class="card tutorial-card">
            <div class="card-header collapsed card-link" data-toggle="collapse" data-target="#collapseThree">
                <div class="d-flex flex-row tutorial-card-header-1">
                    <div class="d-flex flex-row tutorial-card-header-2">
                        <button class="btn btn-dark btn-sm"></button>
                        Interpretable factor model of scVI
                    </div>
                    <span class="badge gs-badge-link">

:doc:`Straight to tutorial...<notebooks/linear_decoder>`

.. raw:: html

                    </span>
                </div>
            </div>
            <div id="collapseThree" class="collapse" data-parent="#accordion">
                <div class="card-body">

It's scVI, but with PCA-like interpretability.


.. image:: ../_static/tutorials/ldvae.svg
   :align: center
   :height: 300px

.. raw:: html

                    <div class="d-flex flex-row">
                        <span class="badge gs-badge-link">

:doc:`To the tutorial <notebooks/linear_decoder>`


.. raw:: html

                        </span>
                    </div>
                </div>
            </div>
        </div>

        <div class="card tutorial-card">
            <div class="card-header collapsed card-link" data-toggle="collapse" data-target="#collapseFour">
                <div class="d-flex flex-row tutorial-card-header-1">
                    <div class="d-flex flex-row tutorial-card-header-2">
                        <button class="btn btn-dark btn-sm"></button>
                        Integration of CITE-seq and scRNA-seq data
                    </div>
                    <span class="badge gs-badge-link">

:doc:`Straight to tutorial...<notebooks/cite_scrna_integration_w_totalVI>`

.. raw:: html

                    </span>
                </div>
            </div>
            <div id="collapseFour" class="collapse" data-parent="#accordion">
                <div class="card-body">

totalVI can be used to integrate datasets from CITE-seq (RNA + protein) and datasets with only RNA (scRNA-seq). Integration enables imputation of missing proteins in the cells measured with scRNA-seq.


.. raw:: html

                    <div class="d-flex flex-row">
                        <span class="badge gs-badge-link">

:doc:`To the tutorial <notebooks/cite_scrna_integration_w_totalVI>`


.. raw:: html

                        </span>
                    </div>
                </div>
            </div>
        </div>

        <div class="card tutorial-card">
            <div class="card-header collapsed card-link" data-toggle="collapse" data-target="#collapseFive">
                <div class="d-flex flex-row tutorial-card-header-1">
                    <div class="d-flex flex-row tutorial-card-header-2">
                        <button class="btn btn-dark btn-sm"></button>
                        Interoperability with R and Seurat
                    </div>
                    <span class="badge gs-badge-link">

:doc:`Straight to tutorial...<notebooks/scvi_in_R>`

.. raw:: html

                    </span>
                </div>
            </div>
            <div id="collapseFive" class="collapse" data-parent="#accordion">
                <div class="card-body">

scvi-tools can be used interfaced directly from R. Learn the basics here.

.. raw:: html

                    <div class="d-flex flex-row">
                        <span class="badge gs-badge-link">

:doc:`To the tutorial <notebooks/scvi_in_R>`

.. raw:: html

                        </span>
                    </div>
                </div>
            </div>
        </div>

        <div class="card tutorial-card">
            <div class="card-header collapsed card-link" data-toggle="collapse" data-target="#collapseScArches">
                <div class="d-flex flex-row tutorial-card-header-1">
                    <div class="d-flex flex-row tutorial-card-header-2">
                        <button class="btn btn-dark btn-sm"></button>
                        Online updates of scvi-tools models via the scArches method
                    </div>
                    <span class="badge gs-badge-link">

:doc:`Straight to tutorial...<notebooks/scarches_scvi_tools>`

.. raw:: html

                    </span>
                </div>
            </div>
            <div id="collapseScArches" class="collapse" data-parent="#accordion">
                <div class="card-body">

scVI, scANVI, and totalVI can be pre-trained on large reference datasets and updated sequentially with query datasets in an online fashion.
This technique uses the method of scArches and enables rapid, and robust integration.

.. raw:: html

                    <div class="d-flex flex-row">
                        <span class="badge gs-badge-link">

:doc:`To the tutorial <notebooks/scarches_scvi_tools>`

.. raw:: html

                        </span>
                    </div>
                </div>
            </div>
        </div>

        <div class="card tutorial-card">
            <div class="card-header collapsed card-link" data-toggle="collapse" data-target="#collapseSix">
                <div class="d-flex flex-row tutorial-card-header-1">
                    <div class="d-flex flex-row tutorial-card-header-2">
                        <button class="btn btn-dark btn-sm"></button>
                        Integration of spatial and sequencing data
                    </div>
                    <span class="badge gs-badge-link">

:doc:`Straight to tutorial...<notebooks/gimvi_tutorial>`

.. raw:: html

                    </span>
                </div>
            </div>
            <div id="collapseSix" class="collapse" data-parent="#accordion">
                <div class="card-body">

gimVI can be used to integrate spatial and sequencing data. Integration enables imputation of missing genes in the cells measured with a spatial technology.


.. raw:: html

                    <div class="d-flex flex-row">
                        <span class="badge gs-badge-link">

:doc:`To the tutorial <notebooks/gimvi_tutorial>`

.. raw:: html

                        </span>
                    </div>
                </div>
            </div>
        </div>

        <div class="card tutorial-card">
            <div class="card-header collapsed card-link" data-toggle="collapse" data-target="#collapseSeven">
                <div class="d-flex flex-row tutorial-card-header-1">
                    <div class="d-flex flex-row tutorial-card-header-2">
                        <button class="btn btn-dark btn-sm"></button>
                        Identification of zero-inflated genes
                    </div>
                    <span class="badge gs-badge-link">

:doc:`Straight to tutorial...<notebooks/AutoZI_tutorial>`

.. raw:: html

                    </span>
                </div>
            </div>
            <div id="collapseSeven" class="collapse" data-parent="#accordion">
                <div class="card-body">

AutoZI can be used to determine which genes are zero-inflated. This can be extended to finding cell-type specific zero-inflation.


.. raw:: html

                    <div class="d-flex flex-row">
                        <span class="badge gs-badge-link">

:doc:`To the tutorial <notebooks/AutoZI_tutorial>`

.. raw:: html

                        </span>
                    </div>
                </div>
            </div>
        </div>

    </div>
    </div>


User-contributed tutorials
--------------------------


.. raw:: html

    <div class="container">
    <div id="accordion" class="shadow tutorial-accordion">

        <div class="card tutorial-card">
            <div class="card-header collapsed card-link" data-toggle="collapse" data-target="#contCollapsedOne">
                <div class="d-flex flex-row tutorial-card-header-1">
                    <div class="d-flex flex-row tutorial-card-header-2">
                        <button class="btn btn-dark btn-sm"></button>
                        Differential expression on Packer C. elegans data
                    </div>
                    <span class="badge gs-badge-link">

:doc:`Straight to tutorial...<notebooks/scVI_DE_worm>`

.. raw:: html

                    </span>
                </div>
            </div>
            <div id="contCollapsedOne" class="collapse" data-parent="#accordion">
                <div class="card-body">

This tutorial was contributed by Eduardo Beltrame.

.. raw:: html

                    <div class="d-flex flex-row">
                        <span class="badge gs-badge-link">

:doc:`To the tutorial <notebooks/scVI_DE_worm>`

.. raw:: html

                        </span>
                    </div>
                </div>
            </div>
        </div>

    </div>
    </div>



.. toctree::
   :maxdepth: 1
   :hidden:

   notebooks/api_overview
   notebooks/data_loading
   notebooks/totalVI
   notebooks/harmonization
   notebooks/linear_decoder
   notebooks/cite_scrna_integration_w_totalVI
   notebooks/scvi_in_R
   notebooks/scarches_scvi_tools
   notebooks/gimvi_tutorial
   notebooks/AutoZI_tutorial
   notebooks/scVI_DE_worm