#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/operators.h>
#include <pybind11/complex.h>
#include <pybind11/functional.h>
#include "Math.h"
#include "Matrix.h"
#include "analysis.h"

namespace py = pybind11;

static RealMatrix RealMatrix_from_nested(const std::vector<std::vector<double>>& a) {
    return RealMatrix(a);
}

static std::vector<std::vector<double>> RealMatrix_to_nested(const RealMatrix& A) {
    std::vector<std::vector<double>> out(A.rows(), std::vector<double>(A.cols()));
    for (size_t i = 0; i < A.rows(); ++i)
        for (size_t j = 0; j < A.cols(); ++j)
            out[i][j] = A.at(i, j);
    return out;
}

void bind_realmatrix(py::module_ &m) {
    py::module_ matrix = m.def_submodule("matrix");

    py::class_<SignedLogDet>(matrix, "SignedLogDet")
        .def_readonly("logdet", &SignedLogDet::logdet)
        .def_readonly("sign", &SignedLogDet::sign)
        .def("__repr__", [](const SignedLogDet& s) {
            return "SignedLogDet(logdet=" + std::to_string(s.logdet) +
                   ", sign=" + std::to_string(s.sign) + ")";
        });

    py::class_<EigenSystem>(matrix, "EigenSystem")
        .def_readonly("D", &EigenSystem::D)
        .def_readonly("P", &EigenSystem::P)
        .def("__repr__", [](const EigenSystem&) { return "EigenSystem(D=..., P=...)"; });


    py::class_<RealMatrix>(matrix, "RealMatrix")
        .def(py::init<>())
        .def(py::init<std::size_t, std::size_t>(), py::arg("rows"), py::arg("cols"))
        .def(py::init<std::vector<double>, std::size_t, std::size_t>(),
             py::arg("data"), py::arg("rows"), py::arg("cols"))
        .def(py::init(&RealMatrix_from_nested), py::arg("data"))

        .def_property_readonly("rows", &RealMatrix::rows)
        .def_property_readonly("cols", &RealMatrix::cols)
        .def_property_readonly("shape", [](const RealMatrix& A) {
            return py::make_tuple(A.rows(), A.cols());
        })

        .def("__getitem__", [](const RealMatrix& A, py::tuple idx) -> double {
            if (idx.size() != 2) throw std::runtime_error("Index must be (i, j)");
            return A.at(idx[0].cast<std::size_t>(), idx[1].cast<std::size_t>());
        })
        .def("__setitem__", [](RealMatrix& A, py::tuple idx, double v) {
            if (idx.size() != 2) throw std::runtime_error("Index must be (i, j)");
            A.at(idx[0].cast<std::size_t>(), idx[1].cast<std::size_t>()) = v;
        })

        .def("at",
             py::overload_cast<size_t, size_t>(&RealMatrix::at),
             py::return_value_policy::reference_internal,
             py::arg("i"), py::arg("j"))
        .def("at",
             py::overload_cast<size_t, size_t>(&RealMatrix::at, py::const_),
             py::arg("i"), py::arg("j"))

        .def("is_symmetric", &RealMatrix::is_symmetric)
        .def("transpose", &RealMatrix::transpose)
        .def("eig", &RealMatrix::eig)
        .def("slogdet", &RealMatrix::slogdet)
        .def("inv", &RealMatrix::inv)

        .def("to_list", &RealMatrix_to_nested)

        .def(-py::self)
        .def(py::self + py::self)
        .def(py::self - py::self)
        .def(py::self * py::self)
        .def(py::self * double())
        .def(double() * py::self)
        .def(py::self / double())

        .def("__repr__", [](const RealMatrix& A) {
            return "RealMatrix(rows=" + std::to_string(A.rows()) +
                   ", cols=" + std::to_string(A.cols()) + ")";
        });

    matrix.def("eye", &eye, py::arg("n"));
    matrix.def("nearest_psd", &nearest_psd, py::arg("R"), py::arg("thr") = 1e-12);
    matrix.def("cholesky_L", &cholesky_L, py::arg("R"));

}

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

    // Math functions
    m.def("Li2", &Li2);
    m.def("Li3", &Li3);
    m.def("CLi2", &CLi2);
    m.def("Cl2", &Cl2);
    m.def("H2", &H2);
    m.def("B", &B);
    m.def("kron", &kron);
    m.def("integrate", &integrate);
    m.def("c_integrate", &c_integrate);

    // Wilson coefficients
    m.def("A0t", &A0t); m.def("F0t", &F0t);
    m.def("B0t", &B0t); m.def("C0t", &C0t);
    m.def("D0t", &D0t); m.def("E0t", &E0t);
    m.def("T", &T);

    m.def("A1t", &A1t); m.def("B1t", &B1t);
    m.def("C1t", &C1t); m.def("D1t", &D1t);
    m.def("E1t", &E1t); m.def("F1t", &F1t);
    m.def("G1t", &G1t);

    m.def("C7t2mt", &C7t2mt); m.def("C7c2MW", &C7c2MW);
    m.def("C8t2mt", &C8t2mt); m.def("C8c2MW", &C8c2MW);
    m.def("F7_1", &F7_1); m.def("F7_2", &F7_2);
    m.def("F8_1", &F8_1); m.def("F8_2", &F8_2);

    m.def("G3H", &G3H); m.def("G4H", &G4H);
    m.def("G7H", &G7H); m.def("G8H", &G8H);
    m.def("EH", &EH);

    m.def("D9H0", &D9H0); m.def("D9H1", &D9H1);
    m.def("Delta3H", &Delta3H); m.def("Delta4H", &Delta4H);
    m.def("Delta7H", &Delta7H); m.def("Delta8H", &Delta8H);
    m.def("C9llH0", &C9llH0); m.def("C9llH1", &C9llH1);
    m.def("C10Wt2mt", &C10Wt2mt); m.def("C10Wc2MW", &C10Wc2MW);
    m.def("C10Zt2mt", &C10Zt2mt); m.def("C10Z2tri", &C10Z2tri);
    m.def("F0SP", &F0SP);

    // Constantes
    m.attr("PI") = PI;
    m.attr("E") = E;
    m.attr("ZETA3") = ZETA3;
    m.attr("GAMMA") = GAMMA;

    // Digamma function
    m.def("psi", &psi);

    // Matrix.h
    py::module matrix = m.def_submodule("matrix");
    py::class_<SparseMatrixWrapper<int>>(matrix, "SparseMatrix")
    .def(py::init<>())
    .def("get_element", &SparseMatrixWrapper<int>::getElement)
    .def("set_element", &SparseMatrixWrapper<int>::setElement)
    .def("get_diagonal_elements", &SparseMatrixWrapper<int>::getDiagonalElements)
    .def("invert", &SparseMatrixWrapper<int>::invert)
    .def("print", &SparseMatrixWrapper<int>::print);

    bind_realmatrix(m);
}
