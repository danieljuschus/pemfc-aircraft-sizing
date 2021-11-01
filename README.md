# Preliminary Propulsion System Sizing Methods for PEM Fuel Cell Aircraft - MSc thesis

This repository contains the main part of the code that I developed for my Master thesis at the TU Delft in 2020/21.

The required Python packages are: ambiance, pyromat, numpy, scipy and cadquery. The latter is only available on Anaconda. Therefore I suggest creating a conda virtual environment with Python 3.8.

The [main file (matlab_engine_sizing_new.py)](matlab_engine_sizing_new.py) is a Python script which sizes a fuel cell 
(FC) system. It is designed to be used by the Initiator, a Matlab  aircraft sizing tool developed at our faculty. 
However, it can also be used directly from Python, for debugging, validation or for sensitivity studies. The main script 
itself calls a number of functions which size various parts of the system. I have also included all the code that I 
used to make all figures in the sensitivity study and validation sections. Furthermore, I have made sure to add a lot of
comments to make everything as clear as possible.

For the methodology applied in the code please take a look at chapter 3 of my thesis report: 
https://repository.tudelft.nl/islandora/object/uuid%3Afdc14875-175e-4a4b-8cde-f2f4b6028194?collection=education. The 
report also contains an introduction into FC aircraft and a summary of the results that were obtained.
