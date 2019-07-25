#include <vector>
#include <iostream>
#include <fstream>

#include "geometry3.hh"
#include "trajectory.hh"
#include "cube.hh"

#include <chrono>

using namespace std;
using namespace std::chrono;

int main (int argc, char **argv) {

  bool errorsfound = false;

  if (argc==1) {
    cout<<"This program needs a command as the first argument.\nList of available commands:\n";
    cout<<"\tsplitgrid: split an integration grid file into the iso- and paratropic parts.\n";
    cout<<"\tgettropplane: get the tropicities of points in a plane perpendicular to either the x, y or z axis.\n";
    cout<<"Run the program without other parameters than the command for operation manual.\n";
    return 7;
  }

  char* command = argv[1];

  if (strcmp(command,"gettropplane") == 0) { //run trop plane getter
//---------------below this line we're checking that the arguments are OK
    if (argc!=7) {
      cout<<"gettropplane expects 5 arguments.\n";
      cout<<"1: Input file in .vti format\n";
      cout<<"2: Magnetic field direction:\n\t0 for x, 1 for -x,\n\t2 for y, 3 for -y,\n\t4 for z, 5 for -z,\n";
      cout<<"3: Fixed coord axis:\n\t0 for x,\n\t1 for y,\n\t2 for z.\n";
      cout<<"4: Fixed coord coordinate as a double. (For example 6.7 or 2.3e2)\n";
      cout<<"5: Output filename\n";
      cout<<"Example: ./main ../QZVPPDh2te-m06-2x.vti 4 2 6.7 output.txt\n";
      return 7;
    }

    try { double dubbel = stod(argv[5]); }
    catch(invalid_argument) { cout<<"The fixed coordinate '"<<argv[5]<<"' was not a double in desired format.\n"; errorsfound = true; }
    
    try { int integer1 = stoi(argv[3]); }
    catch(invalid_argument) { cout<<"The magnetic field direction '"<<argv[3]<<"' was not an integer in desired format.\n"; errorsfound = true; }
    
    try { int integer2 = stoi(argv[4]); }
    catch(invalid_argument) { cout<<"The coord axis '"<<argv[4]<<"' was not an integer in desired format.\n"; errorsfound = true; }
    
    if (stoi(argv[3])<0 || stoi(argv[3])>5){
      cout<<"Expected magnetic field direction:\n\t0 for x, 1 for -x,\n\t2 for y, 3 for -y,\n\t4 for z, 5 for -z,\n";
      cout<<"Entered value was: "<<argv[3]<<"\n";
      errorsfound = true;
    }
    
    if (stoi(argv[4])<0 || stoi(argv[4])>2) {
      cout<<"Expected coord axis:\n\t0 for x,\n\t1 for y,\n\t2 for z.\n";
      cout<<"Entered value was: "<<argv[4]<<"\n";
      errorsfound = true;  
    }
    
    if (errorsfound) {
      cout<<"Run the program with only one parameter - the command name - for operation manual.\n";
      return 7;
    }
//---------------above this line we're checking that the arguments are OK
    
    
    Cube cube(argv[1]);
    vector<vector<int>> tropplane = cube.gettropplane(argv[5], stoi(argv[2]), stoi(argv[3]), stod(argv[4]));
    cube.writetropplane(argv[5],tropplane);
    return 0;
  } else if (strcmp(command,"splitgrid") == 0) {  //run gridsplitter

    if (argc!=6) {
      cout<<"splitgrid expects 4 arguments.\n";
      cout<<"1: Input file in .vti format\n";
      cout<<"2: Magnetic field direction:\n\t0 for x, 1 for -x,\n\t2 for y, 3 for -y,\n\t4 for z, 5 for -z,\n";
      cout<<"3: Gridfile filename\n";
      cout<<"4: Weightfile filename\n";
      cout<<"Example: ./main ../QZVPPDh2te-m06-2x.vti 4 ../grid4h2te ../weights4h2te\n";
      cout<<"The output will be 'gridfile-isotropic', 'weightfile-isotropic', etc...\n";
      return 7;
    } 

    try {int bfielddir=stoi(argv[3]);}
    catch(invalid_argument) { cout<<"The magnetic field direction '"<<argv[3]<<"' couldn't be interpreted as an integer.\n"; return 7;}


    if (stoi(argv[3])<0 || stoi(argv[3])>5){
      cout<<"Expected magnetic field direction:\n\t0 for x, 1 for -x,\n\t2 for y, 3 for -y,\n\t4 for z, 5 for -z,\n";
      cout<<"Entered value was: "<<argv[3]<<"\n";
      return 7;
    }

    Cube cube(argv[2]);
    cube.splitgrid(argv[4],argv[5],stoi(argv[3]));
    return 0;
  } else {
    cout<<"Command not recognized. List of available commands:\n\tgettropplane\n\tsplitgrid\n";
    return 7;
  }

  cout<<"This line was not supposed to be reached in the code.\n";
  return 7;
}
