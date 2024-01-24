import os
from setuptools import setup, find_packages

requires = [
    "pyomo",
    "scipy==1.11",
    "gridx-egret",
    "pint",
]

setup(
    name="gtep",
    version="0.1.dev0",
    python_requires=">=3.7",
    description="a python package",
    packages=find_packages(),
    install_requires=requires,
)
