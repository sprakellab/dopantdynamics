# Dopant Dynamics #

This is the code for the the dopant dynamics project. This project takes place at the Laboratory of Physical Chemistry and Soft
Matter at Wageningen Unversity and Research. The project is part of research performed by Justin Tauber and Ruben Higler 
in the group of Joris Sprakel. 

# Description #

In the simulation project we use HOOMD-BLUE for performing Brownian Dynamics simulations. The code presented here is meant to create/manage input files for these simulations, start simulations via HOOMD-BLUE, data handeling during the simulations and analysis after the simulations have finished.

A new simulation is started by running `hoomd ./dopant.py`. An input file can be generated or given as input.
This will also invoke the bookkeeping process during the simulation (`bookkeeper.py`) and the analysis afterwards (`analyzer.py`).
Basic analysis involves calculating g(r) and mean squared displacements. 

For further analysis functions from the `catpy`-library can be used. Designated modules can be found in the `catpy/ana/`-folder.

To calculate the potential-fields in three dimensions for a perfect crystal and from simulation data several functions are avalaible in `plotpotential.py`.

# dependencies #

* HOOMD-BLUE -- http://glotzerlab.engin.umich.edu/hoomd-blue/

python libraries:
* numpy -- http://www.numpy.org/
* scipy -- http://www.scipy.org/
* scikit-learn -- http://scikit-learn.org/stable/index.html
* mayavi -- http://mayavi.sourceforge.net/
* matplotlib -- http://matplotlib.org/
* lmfit -- https://lmfit.github.io/lmfit-py/

c/c++ programs:
* voro++ -- http://math.lbl.gov/voro++/
* StructureAnalysis for bondparameters -- https://github.com/WolfgangLechner/StructureAnalysis/

Note: to work with the bondparameter functions you should either install StructureAnalysis 
in the same folder as where you run scripts or change the code in the bondparameter module 
to direct to the location of Structure Analysis.Furthermore you will have to change some of 
the code in the StructureAnalysis program with regard to naming and saving files.
