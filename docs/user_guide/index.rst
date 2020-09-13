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
                <img src="../_static/tutorials/totalvi_cell.svg" class="card-img-top" alt="totalvi basic tutorial action icon" height="225">
                <div class="card-body flex-fill">
                    <h5 class="card-title">Analysis of paired RNA and protein data</h5>
                    <p class="card-text">Analyzing CITE-seq data? Use totalVI for joint dimension reduction, integration, differential expression, protein background removal.</p>

.. container:: custom-button

    :doc:`GO<notebooks/totalVI>`


.. raw:: html

                </div>
                </div>
            </div>
            <div class="col-lg-6 col-md-6 col-sm-6 col-xs-12 d-flex">
                <div class="card text-center intro-card shadow">
                <img src="../_static/tutorials/anndata.svg" class="card-img-top" alt="data loading tutorial action icon" height="225">
                <div class="card-body flex-fill">
                    <h5 class="card-title">Data preparation</h5>
                    <p class="card-text">How do I get my data prepared for scvi-tools?</p>

.. container:: custom-button

    :doc:`GO<notebooks/data_loading>`


.. raw:: html

                </div>
                </div>
            </div>
            <div class="col-lg-6 col-md-6 col-sm-6 col-xs-12 d-flex">
                <div class="card text-center intro-card shadow">
                <img src="../_static/tutorials/scanvi.svg" class="card-img-top" alt="scanvi tutorial action icon" height="225">
                <div class="card-body flex-fill">
                    <h5 class="card-title">Cell type annotation</h5>
                    <p class="card-text">Want to transfer cell type annotations to unlabeled cells? scANVI was designed for just that.</p>

.. container:: custom-button

    :doc:`GO<notebooks/harmonization>`


.. raw:: html

                </div>
                </div>
            </div>
            <div class="col-lg-6 col-md-6 col-sm-6 col-xs-12 d-flex">
                <div class="card text-center intro-card shadow">
                <img src="../_static/tutorials/ldvae.svg" class="card-img-top" alt="ldvae tutorial action icon" height="225">
                <div class="card-body flex-fill">
                    <h5 class="card-title">Linearly-decoded scVI</h5>
                    <p class="card-text">It's scVI, but with PCA-like interpretability.</p>

.. container:: custom-button

    :doc:`GO<notebooks/linear_decoder>`


.. raw:: html

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
