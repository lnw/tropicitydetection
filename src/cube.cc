#include <vector>
#include <iostream>
#include <fstream>
#include <list>
#include <iterator>
#include <regex>
#include <cmath>
#include <unistd.h>
// #include <sys/stat.h>

#include "geometry3.hh"
#include "cube.hh"
#include "trajectory.hh"


using namespace std;


// constructor, read from file
Cube::Cube(string filename){
  vector <string> coordinates;
  regex redata("[0-9]E");
  regex reblockstart("<DataArray");
  regex reblockend("</DataArray");
  smatch match;
  bool inblock = 0;
  regex imagedata("<ImageData");
  fstream vti (filename);
  string vtiline;
  if(!vti.good()){
    cout<<"Input file '"<<filename<<"' was not found.\n";
  }

  while (getline (vti, vtiline)) {
    if (inblock == 1 && regex_search(vtiline, match, redata)) {
      coordinates.push_back(vtiline);
      istringstream iss(vtiline);
      vector<string> results((istream_iterator<string>(iss)),istream_iterator<string>());
      coord3d doubleresults(stod(results[0]),stod(results[1]),stod(results[2]));
      field.push_back(doubleresults);
    }
    if (regex_search(vtiline, match, reblockstart)) {
      inblock = 1;
    }
    if (regex_search(vtiline, match, reblockend)) {
      inblock = 0;
    }
    if (regex_search(vtiline, match, imagedata)) {
      istringstream iss(vtiline);
      vector<string> results ((istream_iterator<string>(iss)),istream_iterator<string>());
      xrange = stoi(results[3])+1;
      yrange = stoi(results[5])+1;
      zrange = stoi(results[7])+1;
      origin.push_back(stod(results[10]));
      origin.push_back(stod(results[11]));
      origin.push_back(stod(results[12]));
      spacing.push_back(stod(results[15]));
      spacing.push_back(stod(results[16]));
      spacing.push_back(stod(results[17]));
    }
  }
}


bool Cube::outofbounds (coord3d position) const {
  if (position[0]>xrange || position[1]>yrange || position[2]>zrange || position[0]<0 || position[1]<0 || position[2]<0) {
    return true;
  }
  return false;
}


//linear interpolation
coord3d Cube::getvector(coord3d position) const{
  if (outofbounds(position)) {
    return coord3d(7,7,7); //someone could write a list of these undocumented errorcodes?
  }
  coord3d intpos((int)position[0],(int)position[1],(int)position[2]);
  coord3d sumvec(0,0,0);
  double normsum = 0;
  for (int z=0; z<2; ++z) {
    for (int y=0; y<2; ++y) {
      for (int x=0; x<2; ++x) {
        const double norm = (coord3d(intpos[0]+x,intpos[1]+y,intpos[2]+z)-position).norm();
        if (norm < 1e-12) { // magic number 1e-12: not checking form norm=0 because of floating point inaccuracy
          return field[position[2]*xrange*yrange+position[1]*xrange+position[0]];
        }
        normsum += 1.0/norm;
        sumvec += field[(intpos[2]+z)*xrange*yrange+(intpos[1]+y)*xrange+intpos[0]+x]/norm;
        }
      }
    }

  return sumvec/normsum;
}


//skeleton function for tricubic interpolation later on we figured out that
//it's probably more expensive to use tricubic interpolation than to make up
//for the linear one by increasing grid resolution
coord3d Cube::getvector3(coord3d position) const{
  return coord3d(7,7,7);	
}


vector<vector<TROPICITY>> Cube::gettropplaneZ(double zcoord) const {
  vector<vector<TROPICITY>> tropplaneZ;
  for (int y=0;y<yrange;y++) {
    vector<TROPICITY> point_tropicity;
    tropplaneZ.push_back(point_tropicity);
    for (int x=0;x<xrange;x++){
      trajectory traj(coord3d(x,y,zcoord),getvector(coord3d(x,y,zcoord)),0.01);
      cout<<"\nNEW TRAJECTORY CREATED AT\t"<<x<<","<<y<<","<<zcoord<<"\n";
      traj.complete(*this);
      const string filename = "new-" + to_string(x) + "-" + to_string(y) + "-" + to_string_with_precision(zcoord) + ".txt";
      traj.write2mathematicalist(filename);
      const TROPICITY tr = traj.classify(*this, 4);
      assert(tr != INPUT_ERROR);
      tropplaneZ[y].push_back(tr);
    }
  }
  return tropplaneZ;
}


void Cube::splitgrid(string gridfile, string weightfile, int bfielddir) const{

  vector<coord3d> gridpoints; //coordinates from the grid input file are read into this vector
  vector<double> gridweights; //weights from the weight input file are read into this vector
  vector<string> gridpoints_str; //coordinates from the grid input file are read into this vector
  vector<string> gridweights_str; //weights from the weight input file are read into this vector
  vector<string> dia_points;  //coordinates that were classified as diatropic are written into this vector
  vector<string> dia_weights;  //and corresponding weights into this vector
  vector<string> para_points; //coordinates that were classified as paratropic are written into this vector
  vector<string> para_weights; //and corresponding weights into this vector
  vector<string> zero_points; //if a coordinate couldn't be classified (trajectory got out of bounds), it is written into this vector
  vector<string> zero_weights; //and the corresponding weight into this vector
  vector<string> zero_intensities; // if a coordinate couldn't be classified, the vector at that coord will be written here.
				// lets one check for convergence

  fstream grid (gridfile);
  string gridline;
  if(!grid.good()) { cout<<"Gridfile '"<<gridfile<<"' was not found.\n"; }
  while (getline (grid, gridline)) {
    istringstream gss(gridline);
    vector<string> gridresults((istream_iterator<string>(gss)),istream_iterator<string>());
    coord3d doublegridresults((stod(gridresults[0])-origin[0])/spacing[0],(stod(gridresults[1])-origin[1])/spacing[1],(stod(gridresults[2])-origin[2])/spacing[2]);
    gridpoints.push_back(doublegridresults);
    gridpoints_str.push_back(gridline);
  }

  fstream weights (weightfile);
  string weightsline;
  if(!weights.good()) { cout<<"Weightfile '"<<gridfile<<"' was not found.\n"; }
  while (getline (weights, weightsline)) {
    istringstream wss(weightsline);
    vector<string> weightsresults((istream_iterator<string>(wss)),istream_iterator<string>());
    gridweights.push_back(stod(weightsresults[0]));
    gridweights_str.push_back(weightsline);
  }

//we now have the gridpoints in vector<coord3d> gridpoints and the corresponding weights in vector<double> gridweights
//next: get tropicity at each point, then write out two gridfiles

  for(int i=0;i<gridpoints.size();i++){
    trajectory traj(gridpoints[i],getvector(gridpoints[i]),0.01);
    if (i%100==0){cout<<"i="<<i<<"/"<<gridpoints.size()<<"\n";}
        //cout<<"\nNEW TRAJECTORY CREATED AT\t"<<gridpoints[i]<<"\n";
    traj.complete(*this);
    const TROPICITY classification = traj.classify(*this, bfielddir);
cout << classification << endl;
    assert(classification != INPUT_ERROR);
    if (classification == DIATROPIC){
      dia_points.push_back(gridpoints_str[i]);
      dia_weights.push_back(gridweights_str[i]);
    } else if (classification == PARATROPIC){
      para_points.push_back(gridpoints_str[i]);
      para_weights.push_back(gridweights_str[i]);
    } else if (classification == OUTOFBOUNDS){
      zero_points.push_back(gridpoints_str[i]);
      zero_weights.push_back(gridweights_str[i]);
      ostringstream vectr;
      vectr<<to_string(getvector(gridpoints[i])[0])<<","<<to_string(getvector(gridpoints[i])[2])<<","<<to_string(getvector(gridpoints[i])[2]);
      zero_intensities.push_back(vectr.str());
    } else if (classification == UNCLASSIFYABLE){
      zero_points.push_back(gridpoints_str[i]);
      zero_weights.push_back(gridweights_str[i]);
      ostringstream vectr;
      vectr<<to_string(getvector(gridpoints[i])[0])<<","<<to_string(getvector(gridpoints[i])[1])<<","<<to_string(getvector(gridpoints[i])[2]);
      vectr<<"\t@\t"<<gridpoints_str[i];
      zero_intensities.push_back(vectr.str());
      //cout<<"couldn't classify this point :o(\n";
    }
  }
//now write dia, para and zero points and weights to respective files
  ofstream diapout, diawout, parapout, parawout, zeropout, zerowout, zeroint;
  ostringstream diapoutfile;
  diapoutfile << gridfile << "-diatropic";
  diapout.open(diapoutfile.str());
  for (int i=0;i<dia_points.size();i++) {
    diapout<<dia_points[i]<<"\n";
  }
  diapout.close();
  ostringstream diawoutfile;
  diawoutfile << weightfile << "-diatropic";
  diawout.open(diawoutfile.str());
  for (int i=0;i<dia_weights.size();i++) {
    diawout<<dia_weights[i]<<"\n";
  }
  diawout.close();

  ostringstream parapoutfile;
  parapoutfile << gridfile << "-paratropic";
  parapout.open(parapoutfile.str());
  for (int i=0;i<para_points.size();i++) {
    parapout<<para_points[i]<<"\n";
  }
  parapout.close();
  ostringstream parawoutfile;
  parawoutfile << weightfile << "-paratropic";
  parawout.open(parawoutfile.str());
  for (int i=0;i<para_weights.size();i++) {
    parawout<<para_weights[i]<<"\n";
  }
  parawout.close();

  ostringstream zeropoutfile;
  zeropoutfile << gridfile << "-zerotropic";
  zeropout.open(zeropoutfile.str());
  for (int i=0;i<zero_points.size();i++) {
    zeropout<<zero_points[i]<<"\n";
  }
  zeropout.close();
  ostringstream zerowoutfile;
  zerowoutfile << weightfile << "-zerotropic";
  zerowout.open(zerowoutfile.str());
  for (int i=0;i<zero_weights.size();i++) {
    zerowout<<zero_weights[i]<<"\n";
  }
  zerowout.close();
  zeroint.open("zerointensities.txt");
  for (int i=0;i<zero_intensities.size();i++){
    zeroint<<zero_intensities[i]<<"\n";
  }

}


vector<vector<TROPICITY>> Cube::gettropplane(string filename, int bfielddir, int fixeddir, double fixedcoord) const {
  double steplength = 0.01;
  vector<vector<TROPICITY>> tropplane;
  if (fixeddir==2) {
  fixedcoord = (fixedcoord-origin[1])/spacing[1];
  /// fixedcoord should probably be scaled (according to the .vti header (the gimic outputfile spacing)) at the very first line of this function!
    for (int y=0;y<yrange;y++) {
    cout<<"y = "<<y<<"/"<<yrange<<endl;
    vector<TROPICITY> point_tropicity;
    tropplane.push_back(point_tropicity);
      for (int x=0;x<xrange;x++){
        trajectory traj(coord3d(x,y,fixedcoord),getvector(coord3d(x,y,fixedcoord)),steplength);
        traj.complete(*this);
        const TROPICITY tr = traj.classify(*this, bfielddir);
        assert(tr != INPUT_ERROR);
        tropplane[y].push_back(tr);
      }
    }
  return tropplane;
  }

  else if (fixeddir==1) {
    for (int z=0;z<zrange;z++) {
    cout<<"z = "<<z<<"/"<<zrange<<endl;
    vector<TROPICITY> point_tropicity;
    tropplane.push_back(point_tropicity);
      for (int x=0;x<xrange;x++){
        trajectory traj(coord3d(x,fixedcoord,z),getvector(coord3d(x,fixedcoord,z)),steplength);
        traj.complete(*this);
        const TROPICITY tr = traj.classify(*this, bfielddir);
        assert(tr != INPUT_ERROR);
        tropplane[z].push_back(tr);
      }
    }
  return tropplane;
  }

  else if (fixeddir==0) {
    for (int y=0;y<yrange;y++) {
    cout<<"y = "<<y<<"/"<<yrange<<endl;
    vector<TROPICITY> point_tropicity;
    tropplane.push_back(point_tropicity);
      for (int z=0;z<zrange;z++){
        trajectory traj(coord3d(fixedcoord,y,z),getvector(coord3d(fixedcoord,y,z)),steplength);
        traj.complete(*this);
        const TROPICITY tr = traj.classify(*this, bfielddir);
        assert(tr != INPUT_ERROR);
        tropplane[y].push_back(tr);
      }
    }
  return tropplane;
  }

  else {
    cout<<"FIXEDDIR was not 0-2.\n";
    vector<vector<TROPICITY>> emptyvec;
    return emptyvec;
  }
}


void Cube::writetropplane(string filename, vector<vector<TROPICITY>> tropicities) const{
  ofstream outputfile;
  outputfile.open(filename);
  outputfile<<"trop = {\n";
  for (int i=0;i<tropicities.size();i++){
    outputfile<<tropicities[i];
    if(i<tropicities.size()-1){outputfile<<",";}
    outputfile<<"\n";
  }
  outputfile<<"}";
}


void Cube::writecube(const string& filename) const {
  cout << "File-writing has not yet been implemented. " << filename << "\n";
}


