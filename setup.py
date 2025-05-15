from setuptools import setup, find_packages

setup(
    name="strainscape",
    version="0.1.0",
    packages=find_packages(include=["strainscape", "strainscape.*"]),
    install_requires=[
        "pandas",
        "intervaltree",
        "numpy",
        "snakemake",
    ],
    author="Renee Oles",
    description="A toolkit for analyzing strain evolution in microbiome data",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    python_requires=">=3.6",
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
    ],
) 