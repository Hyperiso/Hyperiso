#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/operators.h>
#include <pybind11/complex.h>
#include "Math.h"


namespace py = pybind11;

void init_math(py::module &m) {

    py::class_<scalar_t>(m, "scalar_t")
        .def(py::init<double, double>(), py::arg("re") = 0.0, py::arg("im") = 0.0)
        .def(py::init<std::complex<double>>())
        .def("real", [](const scalar_t& z) { return z.real(); })
        .def("imag", [](const scalar_t& z) { return z.imag(); })
        .def("to_double", [](const scalar_t& z) { return static_cast<double>(z); })
        .def(py::self + py::self)
        .def(py::self - py::self)
        .def(py::self * py::self)
        .def(py::self / py::self)
        .def(-py::self)
        .def("__repr__", [](const scalar_t& z) {
            return "scalar_t(" + std::to_string(z.real()) + ", " + std::to_string(z.imag()) + ")";
        });


    // Bind scalar_t math functions
    m.def("sqrt", static_cast<scalar_t(*)(const scalar_t&)>(&sqrt));
    m.def("sin", static_cast<scalar_t(*)(const scalar_t&)>(&sin));
    m.def("cos", static_cast<scalar_t(*)(const scalar_t&)>(&cos));
    m.def("tan", static_cast<scalar_t(*)(const scalar_t&)>(&tan));
    m.def("asin", static_cast<scalar_t(*)(const scalar_t&)>(&asin));
    m.def("acos", static_cast<scalar_t(*)(const scalar_t&)>(&acos));
    m.def("atan", static_cast<scalar_t(*)(const scalar_t&)>(&atan));
    m.def("exp", static_cast<scalar_t(*)(const scalar_t&)>(&exp));
    m.def("log", static_cast<scalar_t(*)(const scalar_t&)>(&log));
    m.def("sinh", static_cast<scalar_t(*)(const scalar_t&)>(&sinh));
    m.def("cosh", static_cast<scalar_t(*)(const scalar_t&)>(&cosh));
    m.def("tanh", static_cast<scalar_t(*)(const scalar_t&)>(&tanh));
    m.def("abs", static_cast<scalar_t(*)(const scalar_t&)>(&abs));
    m.def("arg", static_cast<scalar_t(*)(const scalar_t&)>(&arg));
    m.def("norm", static_cast<scalar_t(*)(const scalar_t&)>(&norm));

    // pow overloads
    m.def("pow", static_cast<scalar_t(*)(const scalar_t&, const scalar_t&)>(&pow));
    m.def("pow", static_cast<scalar_t(*)(const scalar_t&, double)>(&pow));
}
