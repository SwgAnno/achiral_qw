# achiral_qw

The aim of this project was to provide a Quantum Walk numerical simulation tool for my B.Sc. thesis project which could support and analyze the behaviour of general Chiral Quantum Walks.
The code was born as a collection of unorganized R scripts, it has now been fully ported to Python and is now shipped as a local package

## Features

The Python project mainly consists in three levels, each one building upon the previous one: the base classes QWGraph and QWGraphCollection, the Simulator classes and analyze methods and lastly all the general purpose plotting and analysis methods.


### Submodules
graph : QWGraph is the basic object for almost any task, it stores and manages all the useful information that are required to simulate a Chiral Quantum Walk on a specific graph

collection : The QWGraphCollection object allow to abstract/multiprocess the analysis routine over an arbitraty set of Graphs. CollectionBuilder create and ships the "most relevant" families of graphs.

simulator : this script contains the family of Schrodinger's equation solver classes needed for Quantum Walk evolution . The basic interface of an abstract Solver is then implemented either relying on QuTip (QutipSESolver) or with eigenvalue decomposition (EigenSESolver)

analyze : The class TransportParameters is a wrapper for the set of figures of merit that we identified as relevant for the analysis of transport. Those parameters can then be employed on a handful of methods that analyze the transport behavour identifying the best events and best phases.

trends , plotter , article : Data extraction and plotting methods over single or relevant families of graph

bessel : Mainly math, nothing really useful eventually

### Examples

The folder "Python/examples" unsurprisingly contains various examples of usage of most of the objects and methods of the package (Work in progress)

## Installation

The project has been bundled into a package put is not available on pip to install, you therefore have to clone the repo and locally install the package.

Navigate to the Python dir and run pip install as follows:

```bash
git clone https://github.com/SwgAnno/achiral_qw
cd achiral_qw/Python
pip install .
```


### Check the version you have installed
You can now check the version of 'achiralqw' you have installed by using the following command:
```bash
pip show achiralqw
```

## Contribution

Although till now this has been a one man project, the code has grown quite big and I'm totally open to contribution.
Just email me at emilio.annoni@studenti.unimi.it

## License(?)

As of 6/1/2023 the whole project is shipped with the DWYLIUP 0.1 License (Do Whatever You Like,It's a University Project).
The Python package you may install has instead a MIT License


