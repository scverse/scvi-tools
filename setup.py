#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup, find_packages

with open("README.rst") as readme_file:
    readme = readme_file.read()

with open("HISTORY.rst") as history_file:
    history = history_file.read()

requirements = [
    "numpy>=1.16.2",
    "torch>=1.0.1",
    "matplotlib>=3.1.2",
    "scikit-learn>=0.20.3",
    "h5py>=2.9.0",
    "pandas>=0.24.2",
    "loompy>=3.0.6",
    "tqdm>=4.31.1",
    "xlrd>=1.2.0",
    "hyperopt==0.1.2",
    "anndata>=0.7",
    "statsmodels",
    "scanpy>=1.4",
    'dataclasses; python_version < "3.7"',  # for `dataclass`
    "scikit-misc",
]

setup_requirements = ["pip>=18.1"]

test_requirements = [
    "pytest>=4.4",
    "pytest-runner>=5.0",
    "flake8>=3.7.7",
    "coverage>=4.5.1",
    "codecov>=2.0.8",
    "black>=19.3b0",
    "nbconvert>=5.4.0",
    "nbformat>=4.4.0",
    "jupyter>=1.0.0",
    "ipython>=7.1.1",
]

extras_requirements = {
    "notebooks": [
        "louvain>=0.6.1",
        "leidenalg>=0.7.0",
        "python-igraph>=0.7.1.post6",
        "colour>=0.1.5",
        "umap-learn>=0.3.10",
        "seaborn>=0.9.0",
        "leidenalg>=0.7.0",
    ],
    "docs": [
        "sphinx>=2.0.1",
        "nbsphinx",
        "sphinx_autodoc_typehints",
        "sphinx-rtd-theme>=0.3.1",
        "autodocsumm",
        "nbsphinx-link",
    ],
    "test": test_requirements,
}
author = (
    "Romain Lopez, "
    "Jeffrey Regier, "
    "Maxime Langevin, "
    "Edouard Mehlman, "
    "Yining Liu, "
    "Achille Nazaret, "
    "Gabriel Misrachi, "
    "Oscar Clivio, "
    "Pierre Boyeau, "
    "Adam Gayoso"
)

setup(
    author=author,
    author_email="romain_lopez@berkeley.edu",
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Natural Language :: English",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: Microsoft :: Windows",
        "Operating System :: POSIX :: Linux",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    description="Single-cell Variational Inference",
    install_requires=requirements,
    license="MIT license",
    long_description=readme + "\n\n" + history,
    include_package_data=True,
    keywords="scvi",
    name="scvi",
    python_requires=">=3.6",
    packages=find_packages(),
    setup_requires=setup_requirements,
    test_suite="tests",
    tests_require=test_requirements,
    extras_require=extras_requirements,
    url="https://github.com/YosefLab/scVI",
    version="0.6.5",
    zip_safe=False,
)
