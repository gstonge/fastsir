/*
 * MIT License
 *
 * Copyright (c) 2022 Guillaume St-Onge
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <BaseContagion.hpp>
#include <DiscreteSIR.hpp>
#include <ContinuousSIR.hpp>

using namespace std;
using namespace fastsir;

namespace py = pybind11;


PYBIND11_MODULE(fastsir, m)
{
    /* ===========
     * Base class
     * ===========*/

    py::class_<BaseContagion>(m, "BaseContagion")

        .def(py::init<EdgeList>(), R"pbdoc(
            Default constructor of the class BaseContagion.

            Args:
               edge_list: Edge list for the network structure.
            )pbdoc", py::arg("edge_list"))

        .def("size", &BaseContagion::size, R"pbdoc(
            Returns the number of nodes.
            )pbdoc")

        .def("get_state_vector", &BaseContagion::get_state_vector, R"pbdoc(
            Returns the vector of state for each node.
            )pbdoc")

        .def("get_current_time", &BaseContagion::get_current_time, R"pbdoc(
            Returns the current time for the process.
            )pbdoc")

        .def("get_number_of_infected_nodes",
                &BaseContagion::get_number_of_infected_nodes, R"pbdoc(
            Returns the number of infected nodes.
            )pbdoc")

        .def("infect_fraction", &BaseContagion::infect_fraction, R"pbdoc(
            Infect a fraction of the nodes.

            Args:
               fraction: Fraction to be infected.
            )pbdoc", py::arg("fraction"))

        .def("infect_node_set", &BaseContagion::infect_node_set, R"pbdoc(
            Infect the nodes in the node set.

            Args:
               node_set: Set of nodes to infect.
            )pbdoc", py::arg("node_set"))

        .def("clear", &BaseContagion::clear, R"pbdoc(
            Recover all nodes.
            )pbdoc")

        .def("reset", &BaseContagion::reset, R"pbdoc(
            Reset time and system, with all susceptible nodes.
            )pbdoc")

        .def("seed", &BaseContagion::seed,
                R"pbdoc(
            Seed the RNG.

            Args:
               seed: seed for the RNG.
            )pbdoc", py::arg("seed"))

        .def("evolve", &BaseContagion::evolve,
                R"pbdoc(
            Let the system evolve over a period of time.

            Args:
               period: Time period of the evolution.
               save_transmission_tree: keep track of transmission
               save_macro_state: keep track of the macro state
            )pbdoc", py::arg("period"), py::arg("save_transmission_tree")=true,
                py::arg("save_macro_state")=true)
        ;


    /* =================================
     * Class deriving from BaseContagion
     * =================================*/

    py::class_<DiscreteSIR, BaseContagion>(m, "DiscreteSIR")

        .def(py::init<EdgeList, double,std::vector<double>>(), R"pbdoc(
            Default constructor of the class DiscreteSIR

            Args:
               edge_list: Edge list for the network structure.
               recovery_probability: Double for the recovery probability
               infection_probability: vector for the infection
                                      probability for different infected
                                      degree.
            )pbdoc", py::arg("edge_list"),
                py::arg("recovery_probability"),
                py::arg("infection_probability"))

        .def("get_lifetime", &DiscreteSIR::get_lifetime, R"pbdoc(
            Returns the lifetime for the current state.
            )pbdoc")

        .def("get_current_macro_state", &DiscreteSIR::get_current_macro_state, R"pbdoc(
            Returns the current macro state
            )pbdoc")

        .def("get_transmission_tree", &DiscreteSIR::get_transmission_tree, R"pbdoc(
            Returns the transmission tree
            )pbdoc")

        .def("get_macro_state_vector", &DiscreteSIR::get_macro_state_vector, R"pbdoc(
            Returns the vector of macro states
            )pbdoc")
        ;

    py::class_<ContinuousSIR, BaseContagion>(m, "ContinuousSIR")

        .def(py::init<EdgeList, double,std::vector<double>>(), R"pbdoc(
            Default constructor of the class ContinuousSIR

            Args:
               edge_list: Edge list for the network structure.
               recovery_rate: Double for the recovery rate
               infection_rate: vector for the infection
                                      rate for different infected
                                      degree.
            )pbdoc", py::arg("edge_list"),
                py::arg("recovery_rate"),
                py::arg("infection_rate"))

        .def("get_lifetime", &ContinuousSIR::get_lifetime, R"pbdoc(
            Returns the lifetime for the current state.
            )pbdoc")

        .def("get_current_macro_state", &ContinuousSIR::get_current_macro_state, R"pbdoc(
            Returns the current macro state
            )pbdoc")

        .def("get_transmission_tree", &ContinuousSIR::get_transmission_tree, R"pbdoc(
            Returns the transmission tree
            )pbdoc")

        .def("get_macro_state_vector", &ContinuousSIR::get_macro_state_vector, R"pbdoc(
            Returns the vector of macro states
            )pbdoc")
        ;


}
