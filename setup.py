#!/usr/bin/env python

from setuptools import setup, Extension
from glob import glob

PKG="pcms"
SRC=PKG+"/src/"
INC=PKG+"/include/"
CMOD="tree_node"

ext_mod = Extension(
    name=PKG+"."+CMOD,
    sources=glob(SRC+"*.c"),
    include_dirs=[INC],
    extra_compile_args=['-g', '-O0']
)

setup(
    name=PKG,
    url="https://github.com/spsvihla/pcms",
    version='1.0',
    description='Phylogenetic Covariance Matrix Sparsification',
    author='Sean Svihla',
    setup_requires=["numpy"],
    install_requires=[
        "numpy",
        "scipy",
        "matplotlib"
    ],
    packages=[PKG],
    ext_modules=[ext_mod],
    include_package_data=True
)
