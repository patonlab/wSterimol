# wSterimol

wSterimol is an automated computational workflow which can be used to obtain multidimensional Sterimol parameters for a conformational ensemble of a given substituent.

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes.

### Prerequisites

**Required**
* [Pymol](https://pymol.org/2/) - Pymol is available for all platforms. Different licenses are available, from commercial to open-source. Open-source pymol is free and contains all the necessary features for wSterimol. However, open-source needs Python as a prerequisite to use it. Several installation guides are available on Internet. It is required as the code uses Pymol API to generate and visualize the conformers.

**MOPAC, Gaussian or ORCA**

* [MOPAC](http://openmopac.net/) - Remember to specify the absolute path in setup.ini for wSterimol to know where to find the software.

* [Gaussian](http://gaussian.com/) - Remember to specify the absolute path in setup.ini for wSterimol to know where to find the software.

* [ORCA](https://orcaforum.cec.mpg.de/) - Remember to specify the absolute path in setup.ini for wSterimol to know where to find the software.

### Example - A tutorial with nButyl

No further installation is required as wSterimol uses the API of Pymol.

We provide herein an example on the published wSterimol to show the reader how it works. For a matter of simplicity, the optimisation will be carried out with Mopac in this example. For the use of the module separately, the reader is invited to read the section about the modules in details.
You need two files to start a calculation: “nbutyl.pdb” and “setup.ini”. These are the coordinates of your molecule and the setup variables for wSterimol.

Open the PDB file with Pymol. A graphical interface should appear. The latter is divided in two: one window corresponds to the command lines, and the other is the 3D representation of your molecule.
The default graphic is not very pretty to visualize small molecules. It can look better by loading "visualize.py” script. To do so, write in the Python/Pymol command line window:

```
run YourPathToTheScript/visualize.py # load script
BallnStick nbutyl # nice ball and sticks
```

The initial structure for the conformation sampling is ready to be used. To start wSterimol, one needs to know the atoms involved in the dihedrals and the primary bond of Sterimol. If one clicks on L (for Label), atom identifiers then ID, the atom numbering should appear. In our example, atom 1 and atom 2 are the primary bond for Sterimol. Then nButyl possesses two dihedrals, each one defined by its 4 atoms.

When the quick analysis is done on yoru side, one can start wSterimol calculation. To do so, one needs to load all the modules in Pymol.

```
run YourPathToTheScript/log.py
run YourPathToTheScript/setup.py
run YourPathToTheScript/ccParse.py
run YourPathToTheScript/generate.py
run YourPathToTheScript/filter_gen.py
run YourPathToTheScript/prepare_file.py
run YourPathToTheScript/optimisation.py
run YourPathToTheScript/filter_opt.py
run YourPathToTheScript/sterimoltools.py
run YourPathToTheScript/sterimol.py
run YourPathToTheScript/weight.py
run YourPathToTheScript/wSterimol.py
```

Then one just needs to write this little command line to get it starts.

```
wSterimol [[id 1, id 2, id 3, id 6],[id 2, id 3, id 6, id 11]], 1, 2
```

The different scripts will be called one after another. If there is an error, it will show up immediately. When the calculation is finished, “wSterimol finished” message should appear. In your working folder, two new files and one folder should have appeared as well. “temp” folder contains all the conformers in PDB format. “weighted.txt” contains the wSterimol values.

One can see how looks like the output from wSterimol calculation. Different values can be retrieved from it. One can get wSterimol values, global energy minimum conformer Sterimol values, range of values within an energy window and their associated Boltzmann weight percentages.

One can also start playing a little bit with the structures by visualizing for instance the Van der Waals surfaces by writing depending on what atomic model you desire:
```
Add_VDW nbutyl, cpk
```
Or
```
Add_VDW nbutyl, bondi
```
One can also visualize all the conformers at the different stages of the calculation by simply using:
```
AddConformers temp
```
Note: "temp” being the default name for the folder containing all the conformers.

## Authors

* **A. V. Brethomé** - *Initial work*
* **S. P. Fletcher** - *Initial work*
* **R. S. Paton** - *Initial work*

## License

This project is licensed under the MIT License.

## Citation

Please, cite this work using:

* To be determined
