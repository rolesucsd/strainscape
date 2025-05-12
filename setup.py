from setuptools import setup, find_packages

setup(
    name="ihmp_analysis",
    version="0.1.0",
    packages=find_packages(),
    install_requires=[
        "pandas",
        "intervaltree",
        "numpy",
    ],
    entry_points={
        'console_scripts': [
            'create-trees=ihmp_analysis.trees:main',
            'summarize-coverage=ihmp_analysis.coverage:main',
        ],
    },
    author="Renee Oles",
    description="Analysis tools for iHMP data",
    python_requires=">=3.6",
) 