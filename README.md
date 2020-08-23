# tropicitydetection

A simple program to determine the tropicity at a set of points in a vector
field.  One possible usage is to investigate magnetically induced electric
currents in molecular systems.


## Usage
Running 'main' without any arguments tells you which commands are available.
Currently there are two commands properly implemented: 'splitgrid' and 'gettropplane'.

## 'splitgrid'
'splitgrid' takes as input a .vti file, a magnetic field direction and a GIMIC gridfile and weightfile.
The program classifies each point of the grid as either diatropic, paratropic or unclassified ("zerotropic").
Example:
```
./main splitgrid QZVPPDh2te-m06-2x.vti 4 grid4h2te weights4h2te
```

## 'gettropplane'
'gettropplane' gets the tropicities of each point in a plane perpendicular to either the x, y, or z axis.
The output can be visualised in Mathematica with the commands
```
'Get["tropplanefile.txt"]'
'ListDensityPlot[trop]'
```

Example: 
```
./main gettropplane QZVPPDh2te-m06-2x.vti 4 2 6.7 output.txt
```


## python interface

There is a thin python interface to most of the functionality.  Compiling it
can be chosen separately.  Two example scripts are located in `bin'; the usage
is the following:

```
python assign_plane.py --input QZVPPDh2te-m06-2x.vti --output output.txt --bfielddir 5 --perpdir 0 --perpcoord 0.0
python splitgrid.py --input QZVPPDh2te-m06-2x.vti --bfielddir 5 --grid grid4h2te --weights weights4h2te
```

## cuda support

The cuda support passes all tests, but is tested only with a device with
compute capability 5.2.  There are three different implementations of the
trajectory completion, two using a double precision vector field and one using
single precision for the vector field (but still doubles for all other floating
point quantities).  In order to compile for other compute capabilities, the
cmake-file needs to be adapted.


## installation

```
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=<DEBUG,RELEASE> -DENABLE_PYTHON=<TRUE,FALSE> -DENABLE_CUDA=<TRUE,FALSE> -DENABLE_OMP=<TRUE,FALSE>
make
```

