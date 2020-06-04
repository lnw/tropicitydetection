#include <string>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "cube.hh"
#include "geometry3.hh"
#include "trajectory.hh"

namespace py = pybind11;

PYBIND11_PLUGIN(libtropicity) {
  py::module mod("shortname", "there could be an explanation here");

  // enum class Tropicity
  py::class_<Tropicity>(mod, "Tropicity");

  // enum class Direction
  py::class_<Direction>(mod, "Direction");

  // class coord3d
  py::class_<coord3d>(mod, "coord3d")
      .def(py::init<float, float, float>());

  // class cube
  py::class_<Cube>(mod, "cube")
      .def(py::init<string>())
      .def("gettropplane", &Cube::gettropplane)     // int, int, double
      .def("writetropplane", &Cube::writetropplane) // string, Plane<TROP>
      .def("writecube", &Cube::writecube)           // string
      .def("outofbounds", &Cube::outofbounds)          //
      .def("splitgrid", &Cube::splitgrid)              // string, string, int
      .def("classify_points", &Cube::classify_points); // vector<coord3d>, Direction

  // class trajectory
  py::class_<Trajectory>(mod, "trajectory")
      .def(py::init<coord3d, coord3d, float>());


  return mod.ptr();
}
