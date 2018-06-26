#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup, find_packages

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

requirements = [
    "numpy>=1.0",
    "torch>=0.4",
    "matplotlib>=2.0",
    "scikit-learn>=0.18",
    "scipy>=1.0",
    "h5py>=2.8",
    "pandas>=0.2",
    "loompy>=2.0",
    "tqdm >= 4"
]

setup_requirements = ['pytest-runner', ]

test_requirements = ['pytest', ]

author = 'Romain Lopez, Jeffrey Regier, Maxime Langevin, Edouard Mehlman, Yining Liu'

setup(
    author=author,
    author_email='romain_lopez@berkeley.edu',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3.6',
        'Operating System :: MacOS :: MacOS X',
        'Operating System :: Microsoft :: Windows',
        'Operating System :: POSIX :: Linux',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
    description="Single-cell Variational Inference",
    install_requires=requirements,
    license="MIT license",
    long_description=readme + '\n\n' + history,
    include_package_data=True,
    keywords='scvi',
    name='scvi',
    packages=find_packages(),
    setup_requires=setup_requirements,
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/YosefLab/scVI',
    version='0.1.3',
    zip_safe=False,
)
