from setuptools import setup

from parnas.version import __version__

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    install_requires=[
        'numpy<1.23,>=1.17',
        'numba==0.55.2',
        'biopython>=1.67',
        'dendropy>=4.5.0',
        'phylo-treetime>=0.9.4'
    ],
    name="parnas",
    version=__version__,
    author="Alexey Markin, Sanket Wagle, Siddhant Grover",
    author_email="alex.markin57@gmail.com",
    license='MIT',
    description="Representative taxon sampling from phylogenetic trees",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/flu-crew/parnas",
    packages=["parnas", "parnas.medoids"],
    classifiers=[
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    entry_points={"console_scripts": ["parnas=parnas.cli:run_parnas_cli"]},
    py_modules=["parnas"],
)
