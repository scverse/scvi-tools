Documentation
********************
``scvi-tools`` (single-cell variational inference tools) [Gayoso21]_ is a package for end-to-end analysis of single-cell omics data primarily developed and maintained by the `Yosef Lab
<https://yoseflab.github.io/>`_ at UC Berkeley. ``scvi-tools`` has two components

* Interface for easy use of a range of probabilistic models for single-cell omics (e.g., scVI, scANVI, totalVI).
* Tools to build new probabilistic models, which are powered by PyTorch, PyTorch Lightning, and Pyro.

If you find a model useful for your research, please consider citing the corresponding publication, which can be found in the corresponding model documentation.

.. panels::
    :card: + intro-card text-center
    :column: col-lg-6 col-md-6 col-sm-6 col-xs-12 d-flex

    ---
    :img-top: _static/computer-24px.svg

    Installation
    ^^^^^^^^^^^^^^^

    New to *scvi-tools*? Check out the installation guide.

    +++

    .. link-button:: installation
            :type: ref
            :text: To the installation guide
            :classes: btn-block btn-secondary stretched-link

    ---
    :img-top: _static/play_circle_outline-24px.svg

    User guide
    ^^^^^^^^^^

    The user guide provides distilled mathematical descriptions of
    the models implemented in scvi-tools and connects the math
    with the code.

    +++

    .. link-button:: user_guide/index
            :type: ref
            :text: To the user guide
            :classes: btn-block btn-secondary stretched-link

    ---
    :img-top: _static/library_books-24px.svg

    API reference
    ^^^^^^^^^^^^^

    The API reference contains a detailed description of
    the scvi-tools API.

    +++

    .. link-button:: api/index
            :type: ref
            :text: To the API reference
            :classes: btn-block btn-secondary stretched-link

    ---
    :img-top: _static/code-24px.svg

    Tutorials
    ^^^^^^^^^^^

    The tutorials walk you through real-world applications of scvi-tools models.
    Developer tutorials help you build new probabilistic models.

    +++

    .. link-button:: tutorials/index
            :type: ref
            :text: To the tutorials
            :classes: btn-block btn-secondary stretched-link


.. toctree::
   :maxdepth: 3
   :titlesonly:
   :hidden:

   installation
   tutorials/index
   user_guide/index
   api/index
   release_notes/index
   references
   contributing/index
   Discussion <https://discourse.scvi-tools.org>
