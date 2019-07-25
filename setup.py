#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup, find_packages

with open("README.rst") as readme_file:
    readme = readme_file.read()

with open("HISTORY.rst") as history_file:
    history = history_file.read()

requirements = [
    "numpy>=1.16.4",
    "torch>=0.4.1",
    "matplotlib>=2.0",
    "scikit-learn>=0.18, <0.20.0",
    "scipy>=1.1",
    "h5py>=2.8",
    "pandas>=0.2",
    "loompy>=2.0",
    "tqdm >= 4",
    "anndata >= 0.6",
    "xlrd >= 1.0",
    "jupyter>=1.0.0",
    "nbconvert>=5.4.0",
    "nbformat>=4.4.0",
    "ipython>=7",
    "umap-learn>=0.3.7",
    "seaborn>=0.9.0",
    "hyperopt>=0.1.2",
]

setup_requirements = ["pytest-runner"]

test_requirements = ["pytest"]

extras_requirements = {
    "test": ["scanpy", "louvain", "leidenalg>=0.7.0", "python-igraph>=0.7.1", "colour"]
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
    packages=find_packages(),
    setup_requires=setup_requirements,
    test_suite="tests",
    tests_require=test_requirements,
    extras_require=extras_requirements,
    url="https://github.com/YosefLab/scVI",
    version="0.4.0",
    zip_safe=False,
)
