[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg?style=flat-square)](https://opensource.org/licenses/MIT)
[![GitHub Workflow Status](https://img.shields.io/github/actions/workflow/status/cda-tum/mqt-misim/ci.yml?branch=main&logo=github&style=flat-square)](https://github.com/cda-tum/mqt-misim/actions?query=workflow%3A%22CI%22)
[![Codecov](https://img.shields.io/codecov/c/github/cda-tum/mqt-misim/main?label=codecov&logo=codecov&style=flat-square)](https://codecov.io/gh/cda-tum/mqt-misim)

<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="https://raw.githubusercontent.com/cda-tum/qmap/main/docs/source/_static/mqt_light.png" width="60%">
    <img src="https://raw.githubusercontent.com/cda-tum/qmap/main/docs/source/_static/mqt_dark.png" width="60%">
  </picture>
</p>

# MQT MiSiM - A Tool for the Simulation of Mixed-Dimensional Quantum Circuits based on Decision Diagrams Written in C++

A tool for the simulation of mixed-dimensional quantum circuit by
the [Chair for Design Automation](https://www.cda.cit.tum.de/) at
the [Technical University of Munich](https://www.tum.de/).

If you have any questions, feel free to contact us via [quantum.cda@xcit.tum.de](mailto:quantum.cda@xcit.tum.de) or by
creating an [issue](https://github.com/cda-tum/mqt-misim/issues) on GitHub. For more information on Decision Diagrams,
please visit [www.cda.cit.tum.de/research/quantum_dd/](https://www.cda.cit.tum.de/research/quantum_dd/).

## Usage

This package caters primarily to our requirements regarding quantum-related functionality and, hence, may not be
straightforward to use for other purposes.

A small example shows how to create set a single qubit in superposition.

```c++
#include "dd/MDDPackage.hpp"

#include <memory>

int main() {
    std::vector<std::size_t> lines{2, 3};
    unsigned int             numLines = 2U;

    auto dd         = std::make_unique<dd::MDDPackage>(numLines, lines); // Create new package instance capable of handling a qubit and a qutrit
    auto zero_state = dd->makeZeroState(numLines);                              // zero_state = |0>

    /* Creating a DD requires the following inputs:
 * 1. A matrix describing a single-qubit/qudit operation (here: the Hadamard matrix)
 * 2. The number of qudits the DD will operate on (here: two lines)
 * 3. The operations are applied to the qubit q0 and the qutrit q1
 * (4. Controlled operations can be created by additionally specifying a list of control qubits before the target declaration)
 */
    auto h_on_qubit  = dd->makeGateDD<dd::GateMatrix>(dd::H(), numLines, 0);
    // auto h_on_qutrit = dd->makeGateDD<dd::TritMatrix>(dd::H3(), 2, 1);

    // Multiplying the operation and the state results in a new state, here a single qubit in superposition
    auto psi = dd->multiply(h_on_qubit, zero_state);

    // Multiplying the operation and the state results in a new state, here a single qutrit in superposition
    // psi = dd->multiply(h_on_qutrit, zero_state);

    // An example of how to create a set of controls and add them together to create a more complex controlled operation
    dd::Controls      control{};
    const dd::Control c{0, 1};
    control.insert(c);

    // An example of a controlled qutrit X operation, controlled on the level 1 of the qubit
    auto CEX = dd->makeGateDD<dd::TritMatrix>(dd::X3, numLines, control, 1);

    psi = dd->multiply(CEX, psi);

    // The last lines retrieves the state vector and prints it
    dd->printVector(psi);
}

```

### System Requirements

Building (and running) is continuously tested under Linux, MacOS, and Windows using
the [latest available system versions for GitHub Actions](https://github.com/actions/virtual-environments). However, the
implementation should be compatible with any current C++ compiler supporting C++17 and a minimum CMake version of 3.19.

### Setup, Build, and Run

To start off, clone this repository using

```shell
git clone --recurse-submodules -j8 https://github.com/cda-tum/mqt-misim/
```

Note the `--recurse-submodules` flag. It is required to also clone all the required submodules. If you happen to forget
passing the flag on your initial clone, you can initialize all the submodules by
executing `git submodule update --init --recursive` in the main project directory.

The MDD (Mixed-Decision-Diagrams) package is a header-only library, so there is no need to build it. However, we provide
a CMake-based build system for testing and benchmarking purposes. Additionally, a dedicated CMake
target `MQT::MiSiM` is provided to easily integrate the DD package into other CMake projects (see,
e.g., https://github.com/cda-tum/qfr).

If you want to build the tests and benchmarks, you need to first configure the project using CMake. This can be done by
calling

```shell
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_MISIM_TESTS=ON
```

This tells CMake to search the current directory `.` (passed via `-S`) for a _CMakeLists.txt_ file and process it into a
directory `build` (passed via `-B`). The flag `-DCMAKE_BUILD_TYPE=Release` tells CMake to configure a _Release_ build (
as opposed to, e.g., a _Debug_ build), while the flag `-DBUILD_MISIM_TESTS=ON` tells CMake to also include the
tests in the build process.

After configuring with CMake, the project can be built by calling

```shell
 cmake --build build --config Release
```

This tries to build the project in the `build` directory (passed via `--build`). Some operating systems and developer
environments explicitly require a configuration to be set, which is why the `--config` flag is also passed to the build
command. The flag `--parallel <NUMBER_OF_THREADS>` may be added to trigger a parallel build.

This generates a number of executables in the `build/test` directory, including the test executable `MiSiM_test`
and the example executable `MiSiM_example`.

## Further Information

The following papers provide further information on different aspects of representing states and operation in the
quantum realm.

- K. Mato, S. Hillmich and R. Wille, "[Mixed-Dimensional Quantum Circuit Simulation with Decision Diagrams](https://www.cda.cit.tum.de/files/eda/2023_mixed_dimensional_quantum_circuit_simulation_with_decision_diagrams.pdf)," 2023 IEEE International Conference on Quantum Computing and Engineering (QCE), Bellevue, Washington, USA, 2023.

- For the representation of unitary matrices and state vectors (with a particular focus on simulation and measurement):
  A. Zulehner and R. Wille. Advanced Simulation of Quantum Computations. IEEE Transactions on Computer Aided Design of
  Integrated Circuits and Systems (TCAD), 2018.
- For the representation and manipulation of unitary matrices (including proof of canonicy, multi-qubit systems, etc):
  P. Niemann, R. Wille, D. M. Miller, M. A. Thornton, and R. Drechsler. QMDDs: Efficient Quantum Function Representation
  and Manipulation. IEEE Transactions on Computer Aided Design of Integrated Circuits and Systems (TCAD), 35(1):86-99, 2016.
- The paper describing this decision diagram package (with a special focus on the representation of complex numbers):
  A. Zulehner, S. Hillmich and R. Wille. How to Efficiently Handle Complex Values? Implementing Decision Diagrams for
  Quantum Computing. The IEEE/ACM International Conference on Computer-Aided Design (ICCAD). 2019
