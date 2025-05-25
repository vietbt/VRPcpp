# VRP C++ Bindings for Python

This repository provides Python C++ bindings for Vehicle Routing Problem (VRP) environments, created using Pybind11. It serves as an extension to the main VRP_াইল research project, enabling the use of highly optimized C++ VRP solvers within a Python environment.

This wrapper facilitates the use of two main VRP solver environments:
* **cvrp_cpp**: For Capacitated Vehicle Routing Problems (CVRP) based on the Hybrid Genetic Search (HGS) algorithm.
* **evrp_cpp**: For Electric Vehicle Routing Problems (EVRP) based on a Variable Neighborhood Search (VNS) approach.

These bindings are part of the research on "Imitation Improvement Learning for Large-Scale Capacitated Vehicle Routing Problems". For more details on the underlying algorithms and research, please refer to the [main project repository](https://github.com/vietbt/VRPpp)  and the [research paper](https://ink.library.smu.edu.sg/sis_research/8025).

## Installation

To install the package, clone this repository and run the following command in the root directory:

```bash
pip install -e .
```

This will build the C++ extensions and install the package in editable mode.

## Overview

The core C++ code implements efficient VRP solvers. This repository wraps that C++ code, exposing it as Python modules (`cvrp_cpp` and `evrp_cpp`) so it can be easily integrated into Python-based research and applications, particularly those involving reinforcement learning or other machine learning approaches for solving VRPs.

The C++ code includes implementations for:
* **CVRP (Capacitated Vehicle Routing Problem)**: Solved using a Hybrid Genetic Search (HGS) algorithm. The `hgs` directory contains the relevant C++ source files (`Genetic.cpp`, `Individual.cpp`, `LocalSearch.cpp`, `Params.cpp`, `Population.cpp`, `Split.cpp`, `main.cpp`, etc.).
* **EVRP (Electric Vehicle Routing Problem)**: Solved using a Variable Neighborhood Search (VNS) based metaheuristic. The `vns` directory contains the `main.cpp` file for this solver.

## Functionality

The Python bindings allow users to:
* Initialize CVRP and EVRP environments.
* Pass problem instances (nodes, capacities, demands, etc.) from Python to the C++ solvers.
* Execute the C++ solvers.
* Retrieve solutions back into the Python environment.

This is particularly useful for frameworks that leverage Python for high-level logic, data manipulation, and machine learning, while relying on C++ for performance-critical computations. The research paper "Imitation Improvement Learning for Large-Scale Capacitated Vehicle Routing Problems" demonstrates such a use case, where these bindings facilitate an imitation learning framework.

## Intended Use

This repository is primarily intended for researchers and developers working on Vehicle Routing Problems who wish to:
* Utilize efficient C++ based VRP solvers within a Python environment.
* Extend or integrate these solvers with machine learning frameworks.
* Reproduce or build upon the research presented in "Imitation Improvement Learning for Large-Scale Capacitated Vehicle Routing Problems".

## Author

* **VietBT** (tvbui@smu.edu.sg)

## Citation

If you use this code or the associated research in your work, please cite:

BUI, The Viet and MAI, Tien. Imitation improvement learning for large-scale capacitated vehicle routing problems. (2023). Proceedings of the 33rd International Conference on Automated Planning and Scheduling (ICAPS 2023): Prague, July 8-13. 1-9. Available at: https://ink.library.smu.edu.sg/sis_research/8025 
