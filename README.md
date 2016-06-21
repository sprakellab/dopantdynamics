# Dopant Dynamics

This folder contains the scripts for the dopant Simulation project.
This project takes place at the Laboratory of Physical Chemistry and Soft
Matter at Wageningen Unversity and Research centre.

The project is part of research performed by Justin Tauber and Ruben Higler 
in the group of Joris Sprakel. 

#dependencies

python libraries:
numpy -- http://www.numpy.org/
scipy -- http://www.scipy.org/
scikit-learn -- http://scikit-learn.org/stable/index.html
mayavi -- http://mayavi.sourceforge.net/
matplotlib -- http://matplotlib.org/
lmfit -- https://lmfit.github.io/lmfit-py/

c/c++ programs:
voro++ -- http://math.lbl.gov/voro++/
StructureAnalysis for bondparameters -- https://github.com/WolfgangLechner/StructureAnalysis/

Note: to work with the bondparameter functions you should either install StructureAnalysis 
in the same folder as where you run scripts or change the code in the bondparameter module 
to direct to the location of Structure Analysis.Furthermore you will have to change some of 
the code in the StructureAnalysis program with regard to naming and saving files.
