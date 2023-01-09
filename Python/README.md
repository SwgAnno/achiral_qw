(github README.md copypaste)

# achiral_qw

The aim of this project was to provide a Quantum Walk numerical simulation tool for my B.Sc. thesis project which could support and analyze the behaviour of general Chiral Quantum Walks.
The code was born as a collection of unorganized R scripts, it has now been fully ported to Python and is ready to be converted into a package

## Features

The project mainly consists in three levels, each one building upon the previous one: the base classes QWGraph and QWGraphCollection, the Simulator and Analyzer classes and lastly all the general purpose plotting and analysis methods.

Graph.py : QWGraph is the basic object for almost any task, it stores and manages all the useful information that are required to simulate a Chiral Quantum Walk on a specific graph

simulator.py : this script contains the Analyzer objec which abstracts over the process of simulating and analyzing transport events on a single graph. Given a graph it can indentify the best transport maximum according to various criteria possibly accounting for all the chiral phases in the process.

collection.py : The QWGraphCollection object allow to abstract/multiprocess the analysis routine over an arbitraty set of Graphs. CollectionBuilder create and ships the "most relevant" families of graphs.

trends.py , plotter.py , article.py : Data extraction and plotting methods over single or relevant families of graph

bessel.py : Mainly math, nothing really useful eventually

### Exapmles

The folder "examples" unsurprizingly contains various examples of usage of all the objects and methods of the package

## Installation

The scripts haven't been bundled into a package yet( though they are planned to be), the repo has to be cloned and the relevan file must be coped into your working directory 

## Contribution

Although till now this has been a one man project, the code has grown quite big and I'm totally open to contribution.
Just email me at emilio.annoni@studenti.unimi.it

## License(???)

As of 6/1/2023 the project is shipped with the DWYLIUP 0.1 License (Do Whatever You Like,It's a University Project)


