## **Multi-scale Anisotropic Rough Surface** (MARS) Algorithm 

Thomas O. Jelly / 2018.

This package provides a Multi-scale Anisotropic Rough Surface algorithm for numerical generation of rough surfaces with specified statistical parameters. It comes with both CLI and GUI.

## Getting started

### Prerequisites

The algorithm has been written to work only with Python 2. It was developed in Python 2.7.15 environment but any 2.7.X version should run it correctly. Packages used in the project have been listed in ```requirement.txt``` 
Running 
```
pip install -r requirements.txt
``` 
should automatically fetch and install all the necesseray packages. If you are using dual Python environment, then you should be careful to install those packages for the Python2 interpreter by running:
```
python2 -m pip install -r requirements.txt
```

### Running MARS

The program comes with both command-line and graphical interfaces. To run the CLI, just run:
```
python main.py
```
All the parameters of the surface generation are stored in the file. If you want to change any of them, you should edit the file manually as argument parsing has not been added yet.

You can run the GUI by executing the GUI.py file:
```
python GUI.py
```
It provides more detailed control of the surface generation process, as well as more advanced options to the user.

### UNDER CONSTRUCTION!
