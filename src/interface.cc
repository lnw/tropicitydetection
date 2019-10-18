#include <string>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "geometry3.hh"
#include "cube.hh"
#include "trajectory.hh"

namespace py = pybind11;

PYBIND11_PLUGIN(libtropicity) {
  py::module m("lalala", "pybind11 ... whatever");

// class coord3d
  py::class_<coord3d>(m, "coord3d")
    .def(py::init<float, float, float>());

// class cube
  py::class_<Cube>(m, "cube")
    .def(py::init<string>())
//    .def("gettropplaneZ", &Cube::gettropplaneZ) // double
    .def("gettropplane", &Cube::gettropplane) // double
//    .def("writetropplaneZ", &Cube::writetropplaneZ) // string, vector<vector<int>>
    .def("writecube", &Cube::writecube) // string
    .def("getvector", &Cube::getvector) // coord3d
    .def("outofbounds", &Cube::outofbounds) //coord3d
    .def("splitgrid", &Cube::splitgrid); //coord3d

// class trajectory
  py::class_<trajectory>(m, "trajectory")
    .def(py::init<coord3d, coord3d, float>());



    return m.ptr();
}

