# achiral_qw

The aim of this project was to provide a Quantum Walk numerical simulation tool for my B.Sc. thesis project which could support and analyze the behaviour of general Chiral Quantum Walks.
The code was born as a collection of unorganized R scripts, it has now been fully ported to Python and is now shipped as a local package

## Features

The Python project mainly consists in three levels, each one building upon the previous one: the base classes QWGraph and QWGraphCollection, the Simulator and Analyzer classes and lastly all the general purpose plotting and analysis methods.


### Submodules
graph : QWGraph is the basic object for almost any task, it stores and manages all the useful information that are required to simulate a Chiral Quantum Walk on a specific graph

collection : The QWGraphCollection object allow to abstract/multiprocess the analysis routine over an arbitraty set of Graphs. CollectionBuilder create and ships the "most relevant" families of graphs.

simulator : this script contains the Analyzer objec which abstracts over the process of simulating and analyzing transport events on a single graph. Given a graph it can indentify the best transport maximum according to various criteria possibly accounting for all the chiral phases in the process.

trends , plotter , article : Data extraction and plotting methods over single or relevant families of graph

bessel : Mainly math, nothing really useful eventually

### Examples

The folder "Python/examples" unsurprisingly contains various examples of usage of mpst of the objects and methods of the package

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

### Contribution

Although till now this has been a one man project, the code has grown quite big and I'm totally open to contribution.
Just email me at emilio.annoni@studenti.unimi.it

## License(?)

As of 6/1/2023 the whole project is shipped with the DWYLIUP 0.1 License (Do Whatever You Like,It's a University Project).
The Python package you may install has instead a MIT License


