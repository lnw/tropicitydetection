# tropicitydetection
Splitting .vti files to paratropic and diatropic parts.
# Usage
Running 'main' without any arguments tells you which commands are available.
Currently there are two commands properly implemented: 'splitgrid' and 'gettropplane'.
# 'splitgrid'
'splitgrid' takes as input a .vti file, a magnetic field direction and a GIMIC gridfile and weightfile.
The program classifies each point of the grid as either diatropic, paratropic or unclassified ("zerotropic").
Example: './main splitgrid ../QZVPPDh2te-m06-2x.vti 4 ../grid4h2te ../weights4h2te'
# 'gettropplane'
'gettropplane' gets the tropicities of each point in a plane perpendicular to either the x, y or z axis.
The output can be visualized in Mathematica with the commands
'Get["tropplanefile.txt"]'
'ListDensityPlot[trop]'
Example: ./main gettropplane ../QZVPPDh2te-m06-2x.vti 4 2 6.7 output.txt
