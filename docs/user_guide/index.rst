User guide
==========

|Open In Colab|

.. |Open In Colab| image:: https://colab.research.google.com/assets/colab-badge.svg
    :target: https://colab.research.google.com/github/yoseflab/scVI/blob/stable


The easiest way to get familiar with scvi-tools is to follow along with our tutorials.
The tutorials are accessible on the sidebar to the left. They are also designed to work seamlessly in Google Colab, a free cloud computing platform. These tutorials have a Colab badge in their introduction. In general, these tutorials are designed to work with the latest installable version of scvi-tools.

To download the tutorials:

1. Click the Colab within the tutorial (if available).
2. Download it with the option in the file menu.


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
                        scRNA-seq atlas-level integration and label transfer with scVI and scANVI
                    </div>
                    <span class="badge gs-badge-link">

:doc:`Straight to tutorial...<notebooks/harmonization>`

.. raw:: html

                    </span>
                </div>
            </div>
            <div id="collapseTwo" class="collapse" data-parent="#accordion">
                <div class="card-body">

TODO

.. image:: ../_static/tutorials/scanvi.svg
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

    </div>
    </div>


User-contributed tutorials
--------------------------


.. raw:: html

    <div class="container">
    <div id="accordion" class="shadow tutorial-accordion">

        <div class="card tutorial-card">
            <div class="card-header collapsed card-link" data-toggle="collapse" data-target="#collapseOne">
                <div class="d-flex flex-row tutorial-card-header-1">
                    <div class="d-flex flex-row tutorial-card-header-2">
                        <button class="btn btn-dark btn-sm"></button>
                        Differential expression on Packer C. elegans data
                    </div>
                    <span class="badge gs-badge-link">

:doc:`Straight to tutorial...<notebooks/contributed/scVI_DE_worm>`

.. raw:: html

                    </span>
                </div>
            </div>
            <div id="collapseOne" class="collapse" data-parent="#accordion">
                <div class="card-body">

.. raw:: html

                    <div class="d-flex flex-row">
                        <span class="badge gs-badge-link">

:doc:`To the tutorial <notebooks/contributed/scVI_DE_worm>`

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
   notebooks/cite_scrna_integration_w_totalVI
   notebooks/harmonization
   notebooks/AutoZI_tutorial
   notebooks/gimvi_tutorial
   notebooks/linear_decoder
   autotune
   notebooks/contributed/scVI_DE_worm
