# HPE-FEM-PLAYGROUND

## Introduction

This project is a playground on the basis of the famous open-source FEM lib MFEM (`https://github.com/mfem/mfem.git`), and serve the following founctionalities:

- Convert COMSOL mesh (.mphtxt, ./mphbin(TODO)) to standard gmsh file (.msh) or MFEM custon mesh format (.mesh)
- Make abstracts, adapters, decorators and plugins for MFEM, so we may plugin our custom algorithm or data schema for each steps (Assembly, BC Application, Linear solving...(Although I don't think it is easy to even get closed to the performance of the modern and highly optimized Linear Solvers)).
- We'd NEVER touch the following from MFEM:

- - The Mesh data container, and representation of weak form terms / boundary conditions.
- - The vector/matrix and their computation logics. I don't think that I'd like to do better hardware level gemm acceleration.

- We WILL look forward to do:

- - JIT for certain types of problems: The MFEM is for runtime generality, yet for a fixed type of problem we may pre-generate and compile specific C code for better runtime performance.
- - Data locality optimization: from traditional Array Of Objects to the per-data-per-array layout, for better data locality and better Cache hit rate.
- - Combining with ML: For the extreme on-the-fly performance for some type of problems, we may train MPL or even literally and AI to predict the connections between given conditions and physics field values.

## File Structure

``` bash
├── README.md
├── CMakeLists.txt
├── benchmark
├── example
├── results
├── scripts
├── share
├── src
└── testdata
```

- `CMakelists.txt`: Major cmake script.
- `benchmark`: To do benchmarks between my custom plugins with the MFEM's original performance.
- `example`: For my learning of basic use of MFEM (maybe to provide examples for my custom adapters in the future).
- `scripts`: Some small Python scriptsfor some common utilities like mesh format conversion.
- `share`: The cmake configuration will automatically clone all the needed third-party libs to the dir. Use shouldn't manually edit it.
- `src`: To place my custom designs.
- `testdata`: Data for example and benchmark use.
- `results`: To store outputs. Without special commandline args, program results would automatically be stored in subdir of `results`

## How To Run

We highly recommand that you use gcc/g++ version >=11, or you may face the risk of failing to compile the mfem lib.
Also we recommand to divide debug build and release build (Debug build for problem checking and validation, release build for benchmark and release). For example:

``` bash
# To Build for debug
cmake -S . -B cmake-build-debug -DCMAKE_BUILD_TYPE=Debug
cmake --build -j4 cmake-build-debug
# To Build for release
cmake -S . -B cmake-build-release -DCMAKE_BUILD_TYPE=Release
cmake --build -j4 cmake-build-release
```

All dependencies will be automatically managed with the cmake script, so feel free to run the project anywhere.
