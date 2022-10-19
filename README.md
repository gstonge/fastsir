# fastsir

Fast c++ implementation of simulations for nonlinear SIR model on networks

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

See the demos on phase transitions and transmission trees.
