# python-fds-tools
A Python Package to Assist in Developing and Post-Processing Data Produced Through the Computational Fluid Dynamics Software Fire Dynamics Simulator.

This software is part of a python library to assist in developing and analyzing simulation results from Fire Dynamics Simulator (FDS). FDS is an open source computational fluid dynamics (CFD) software package developed by NIST. The sourcecode is available at: https://github.com/firemodels/fds

Installation

This module was developed for use in an anaconda virtual environment. The installation steps listed below use this approach.

*The first step is to clone this repository to your local machine:

git clone https://github.com/johodges/pyfdstools

*The next step is to set up the anaconda virtual environment:

conda env update --file pyfdstools/env.yaml
conda activate fdsTools

*Add the new repository to the python path. This can be done by updating the user path, or adding environmental variables. Note, if this is done through anaconda, the command below needs to be run from the directory from which the repository was cloned, not inside the pyfdstools directory.

conda develop .

*Run the example cases

python pyfdstools/examples.py
