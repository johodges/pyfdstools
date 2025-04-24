# python-fds-tools
A Python Package to Assist in Developing and Post-Processing Data Produced Through the Computational Fluid Dynamics Software Fire Dynamics Simulator.

This software is part of a python library to assist in developing and analyzing simulation results from Fire Dynamics Simulator (FDS). FDS is an open source computational fluid dynamics (CFD) software package developed by NIST. The sourcecode is available at: https://github.com/firemodels/fds

# Installation

This module was developed for use in a virtual environment. 
The package can be installed through pip or source.

* Configuring the virtual environment
  ```
  python -m venv c:\path\to\myenv
  ```
* Activate the virtual environment
  ```
  source c:\path\to\myenv\Scripts\activate
  ```
* (Option 1) Installing via pip
  - Install
    ```
    python -m pip install pyfdstools
    ```
  - Run the example cases
    ```
    python -c 'import pyfdstools as fds; fds.runExamples()'
    ```
* (Option 2) Installing from source
  - Navigate to desired installation location
  - Clone the repository
    ```
    git clone https://github.com/johodges/pyfdstools
    ```
  - Install with pip
    ```
    pip install pyfdstools/
    ```
  - Run the example cases
    ```
    cd pyfdstools
    python pyfdstools/examples.py
    ```
* (Option 3) Installing from source in edit mode. The advantage of this mode is python will import modules from the source directly. This allows you to modify scripts as needed or update with a git pull.
  - Navigate to desired installation location
  - Clone the repository
    ```
    git clone https://github.com/johodges/pyfdstools
    ```
  - Install with pip
    ```
    pip install -e pyfdstools/
    ```
  - Run the example cases
    ```
    cd pyfdstools
    python pyfdstools/examples.py
    ```

<!---
# Integration with BlenderFDS

This module was developed for use in an anaconda virtual environment. The installation steps listed below use this approach.

* Create the anaconda environment

```
conda create -n blenderfds python=3.7.0
```

* Activate the anaconda environment

```
conda activate blenderfds
```

* Update the anaconda environment

```
conda env update --file pyfdstools/blenderfds_env.yaml
```

* Add the new repository to the python path. This can be done by updating the user path, or adding environmental variables. Note, if this is done through anaconda, the command below needs to be run from the directory from which the repository was cloned, not inside the pyfdstools directory.

```
conda develop .
```

* Set environmental variable for QT

```
conda env config vars set QT_PLUGIN_PATH="$CONDA_PREFIX/Library/plugins"
```

* Set environmental variable for Blender installation

```
export BLENDER_INSTALLATION_DIRECTORY="/c/Program\ Files/Blender\ Foundation/Blender\ 2.91/2.91"
```

* Remove pre-installed blender python

```
mv "$BLENDER_INSTALLATION_DIRECTORY/python/" "$BLENDER_INSTALLATION_DIRECTION/_python/"
mv /c/Program\ Files/Blender\ Foundation/Blender\ 2.91/2.91/python/ /c/Program\ Files/Blender\ Foundation/Blender\ 2.91/2.91/_python/
```

* Create symbolic link to anaconda environment

```
ln -s $CONDA_PREFIX "$BLENDER_INSTALLATION_DIRECTORY/python"
ln -s /c/ProgramData/Anaconda3/envs/blenderfds /c/Program\ Files/Blender\ Foundation/Blender\ 2.91/2.91/python
```

* Add blender to the path
* Restart command line
* Activate the blenderfds environment
* Start blender from the command line
--->

# Usage example

This simple example demonstrates how to use pyfdstools to work with 2-D slice data. The example does the following:
* Reads 2-D slice data from the examples directory
* Dumps the output to a csv file in the current working directory
* Visualizes the queried slice in a figure
* Saves the figure to the current working directory

```python
import pyfdstools as fds
import os

# configure case information
# sets the working directory to be the examples directory provided with pyfdstools.
# Note that the working_dir is set here to a zip file
working_dir = os.path.join(os.path.dirname(fds.__file__),'examples','data','case001.zip') 
chid = 'case001' # sets the fds CHID for the case
quantity = 'TEMPERATURE' # sets the quantity to be queried
axis, value = 1, 2.55 # sets the axis (1 = x, 2 = y, 3 = z) and coordinate of the 2-D slice to query
time, dt = 30, 60 # sets the query time and window
qnty_mn, qnty_mx = 0, 1000
cbarnumticks = 11
outdir = os.getcwd() # sets the output directory to the current working directory

# Read 2-D slice data
data, unit = fds.query2dAxisValue(working_dir, chid, quantity, axis, value, time=time, dt=dt)

# Dump 2-D slice data to a csv file
fds.renderSliceCsvs(data, chid, outdir)

# Visualize 2-D slice data and save to a file
fig, ax = fds.plotSlice(data['x'], data['z'], data['datas'][:, :, -1], axis, clabel="%s (%s)"%(quantity, unit), qnty_mn=qnty_mn, qnty_mx=qnty_mx, cbarnumticks=cbarnumticks)

fig.savefig(os.path.join(outdir, '%s_%s_%0.0f_%0.4f_final_frame.png'%(chid, quantity, axis, value)))
fig.show()
```

# Additional examples

The examples directory has a series of example files showing how a user can use pyfdstools to interact with FDS input files and results. All the examples can be run using the runExamples routine:

```python
import pyfdstools as fds
fds.runExamples()
```

The outputs from these examples will be located in the pyfdstools/generated directory.

The examples can also be run individually by running the individual files in the examples directory. The examples can be modified to work with different FDS results by editing the scripts, or by setting some of the parameters with command line arguments. The example below changes the dump_2d_slice_to_csv example to work with case002 instead of case001.

```python
python dump_2d_slice_to_csv.py --chid case002 --quantity TEMPERATURE --axis 3 --value 7.2 --time 30 --dt -1 --working_dir data/case002.zip
```

# Citation

If you use this software in your research, please consider citing this project as:

Hodges, J. L., pyFDStools: A Python Package to Assist in Developing and Post-Processing Data Produced Through the Computational Fluid Dynamics Software Fire Dynamics Simulator, (2020), GitHub repository, https://github.com/johodges/pyfdstools.
