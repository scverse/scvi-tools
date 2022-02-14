Documentation
********************
``scvi-tools`` (single-cell variational inference tools) [Gayoso22]_ is a package for end-to-end analysis of single-cell omics data primarily developed and maintained by the `Yosef Lab
<https://yoseflab.github.io/>`_ at UC Berkeley. ``scvi-tools`` has two components

* Interface for easy use of a range of probabilistic models for single-cell omics (e.g., scVI, scANVI, totalVI).
* Tools to build new probabilistic models, which are powered by PyTorch, PyTorch Lightning, and Pyro.

If you find a model useful for your research, please consider citing the scvi-tools manuscript as well as the corresponding publication, which can be found in the corresponding model documentation.

.. card:: Installation :octicon:`plug;1em;`
    :link: installation
    :link-type: doc

    New to *scvi-tools*? Check out the installation guide.

.. card:: User guide :octicon:`play;1em;`
    :link: user_guide/index
    :link-type: doc

    The user guide provides distilled mathematical descriptions of
    the models implemented in scvi-tools and connects the math
    with the code.

.. card:: API reference :octicon:`book;1em;`
    :link: api/index
    :link-type: doc

    The API reference contains a detailed description of
    the scvi-tools API.

.. card:: Tutorials :octicon:`workflow;1em;`
    :link: tutorials/index
    :link-type: doc

    The tutorials walk you through real-world applications of scvi-tools models.
    Developer tutorials help you build new probabilistic models.

.. card:: Discussion :octicon:`megaphone;1em;`
    :link: https://discourse.scvi-tools.org

    Need help? Reach out on our forum to get your questions answered!


.. card:: GitHub :octicon:`mark-github;1em;`
    :link: https://github.com/yoseflab/scvi-tools

    Find a bug? Interested in improving scvi-tools? Checkout our GitHub for the latest developments.


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
   GitHub <https://github.com/YosefLab/scvi-tools>
