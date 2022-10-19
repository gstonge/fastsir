# fastsir

Fast c++ implementation of simulations for nonlinear SIR model on networks

## Requirements and dependencies

* A compiler with C++17 support
* `python3`
* if not already installed, `pybind11` and `setuptool` will be installed in an isolated environment
  to build the package.

## Installation

First, clone this repository. If using SSH (recommended)
```bash
git clone git@github.com:gstonge/fastsir.git
```
or

```bash
git clone https://github.com/gstonge/fastsir.git
```

Second, use pip to install the module.
```bash
pip install ./SamplableSet
```
It is recommended to perform the installation in a virtual environment.

## A peak under the hood

On the C++ side, we have a hierarchy of classes inheriting from the base (dummy) class
```
├── src
    ├── BaseContagion.hpp
```

Right now, there is only one derived class for discrete-time SIR models
```
├── src
    ├── DiscreteSIR.hpp
```

The C++ classes are "exposed" to python using [pybind11](https://pybind11.readthedocs.io/en/stable/index.html).
```
├── src
    ├── bind.hpp
```

One should not have to deal with the C++ side, just import the class in python using
```
from fastsir import DiscreteSIR as sir
```

## Basic usage

Everything revolves around the `sir` object (shorthand here for `DiscreteSIR`).
To initialize the object, we need an edge list, a recovery probability for infected nodes, and a
vector of infection probability. The latter gives the probability of getting infected per time step
as a function of the infected degree, (number of infected nodes in the neighborhood).

Let's use the networkx package to get a graph

```
from fastsir import DiscreteSIR as sir
import numpy as np
import networkx as nx

#Network parameters
N = 10000
m = 10
p = 0.2
G = nx.watts_strogatz_graph(N,m,p)
max_degree = np.max([G.degree(n) for n in G])
edgelist = list(G.edges())
```

For the contagion process, we use some arbitrary parameter. For the infection probability, we assume
a "simple contagion"
```
#contagion parameters
r = 0.01 #recovery probability
q = 0.0015 #from the demo on phase transition, we know this is close to criticality
recovery_probability = r
infection_func = lambda k: 1-(1-q)**k
infection_probability = infection_func(np.arange(max_degree+1))
```

Finally, or stochastic process is initialized like this
```
#initialize object
process = sir(edgelist,recovery_probability,infection_probability)
```

Now we need to infect a certain number of nodes initially at random. Here we use a single node, but
more could be used. It must be a set.
```
initial_infected_nodes = {np.random.randint(N)}
process.infect_node_set(initial_infected_nodes)
```

Now to advance in time, we specify the `period` for which we want the process to evolve and call the
method `evolve`:
```
period = 10
process.evolve(period,save_transmission_tree=True,save_macro_state=True)
```
Note: here we can specify if we want to keep track of the transmission tree and of the macro state
(number of susceptible,infected,recovered) at each time.
If we want to evolve until there are no more infected nodes---the process has died---, simply set
`period = np.inf`.

After that, we can probe our process using various methods:
```
#get the current state of the network
S,I,R = process.current_macro_state()

#get the vector of macro state at each time step
macro_state_vector = process.get_macro_state_vector()

#get the transmission tree
tree = process.get_transmission_tree()
```

For more advanced usage, see the demos on transmission trees and phase transitions.
Also, for a complete list of useful methods, do
```
help(sir)
```
or
```
help(process)
```
