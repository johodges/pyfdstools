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
  ```
  python -m pip install pyfdstools
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
* Run the example cases
```
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

# Citation

If you use this software in your research, please consider citing this project as:

Hodges, J. L., pyFDStools: A Python Package to Assist in Developing and Post-Processing Data Produced Through the Computational Fluid Dynamics Software Fire Dynamics Simulator, (2020), GitHub repository, https://github.com/johodges/pyfdstools.
