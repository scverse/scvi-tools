User guide
==========

|Open In Colab|

.. |Open In Colab| image:: https://colab.research.google.com/assets/colab-badge.svg
    :target: https://colab.research.google.com/github/yoseflab/scVI/blob/stable


The easiest way to get familiar with scVI is to follow along with our tutorials!
The tutorials are accessible on the sidebar to the left. Some are designed to work seamlessly in Google Colab, a free cloud computing platform. These tutorials have a Colab badge in their introduction. In general, these tutorials are designed to work with the latest installable version of scVI. Previous versions of the tutorials are available by changing the Read the Docs version (available at the bottom of the page).

To download the tutorials:

1. Click the Colab badge above
2. Open the tutorial
3. Download it with the option in the file menu.
4. When you execute the notebook yourself, please set your own `save_path`.

.. raw:: html

    <div class="container">
        <div class="row">
            <div class="col-lg-6 col-md-6 col-sm-6 col-xs-12 d-flex">
                <div class="card text-center intro-card shadow">
                <img src="../_static/tutorials/scvi_batch.png" class="card-img-top" alt="scvi basic tutorial action icon" height="225">
                <div class="card-body flex-fill">
                    <h5 class="card-title">Introduction to scVI</h5>
                    <p class="card-text">Dimension reduction, integration, differential expression.
                    </p>

.. container:: custom-button

    :doc:`GO<basic_tutorial>`


.. raw:: html

                </div>
                </div>
            </div>
            <div class="col-lg-6 col-md-6 col-sm-6 col-xs-12 d-flex">
                <div class="card text-center intro-card shadow">
                <img src="../_static/tutorials/totalvi_cell.svg" class="card-img-top" alt="totalvi basic tutorial action icon" height="225">
                <div class="card-body flex-fill">
                    <h5 class="card-title">Introduction to totalVI</h5>
                    <p class="card-text">Analyzing CITE-seq data? Joint dimension reduction, integration, differential expression, protein background removal.</p>

.. container:: custom-button

    :doc:`GO<totalvi>`


.. raw:: html

                </div>
                </div>
            </div>
            <div class="col-lg-6 col-md-6 col-sm-6 col-xs-12 d-flex">
                <div class="card text-center intro-card shadow">
                <img src="../_static/tutorials/anndata.svg" class="card-img-top" alt="data loading tutorial action icon" height="225">
                <div class="card-body flex-fill">
                    <h5 class="card-title">Data preparation</h5>
                    <p class="card-text">How do I get my data prepared for scvi?</p>

.. container:: custom-button

    :doc:`GO<data_loading>`


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

    :doc:`GO<harmonization>`


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

    :doc:`GO<linear_decoder>`


.. raw:: html

                </div>
                </div>
            </div>
        </div>
    </div>

.. toctree::
   :maxdepth: 1
   :hidden:

   notebooks/basic_tutorial
   notebooks/data_loading
   notebooks/totalVI
   cite_scrna_integration_w_totalVI
   notebooks/harmonization
   notebooks/AutoZI_tutorial
   notebooks/gimvi_tutorial
   notebooks/linear_decoder
   autotune
   scVI_DE_worm
