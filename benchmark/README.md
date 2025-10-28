# Benchmark Suite for HPC-FEM

This directory contains benchmarking tools for evaluating solver performance and scalability.

## Contents

- `poisson_scaling/` - Poisson equation benchmark application
- `run_scaling.py` - Automation script for scaling studies
- `plot_results.py` - Visualization script for results

## Quick Start

### 1. Build the Benchmark Application

```bash
cd /path/to/hpc-fem-playground
mkdir -p cmake-build-release
cd cmake-build-release
cmake -DCMAKE_BUILD_TYPE=Release ..
make benchmark_poisson -j$(nproc)
```

### 2. Run a Single Benchmark

```bash
cd cmake-build-release/benchmark/poisson_scaling

# Serial AMG solver
./benchmark_poisson --solver amg --mesh-size 32 --order 2

# Parallel DDM solver
mpirun -np 4 ./benchmark_poisson --solver ddm --mesh-size 64 --order 2
```

### 3. Run Scaling Studies

#### Strong Scaling (Fixed Problem Size)

```bash
cd benchmark
python3 run_scaling.py \
    --study strong \
    --mesh-size 64 \
    --order 2 \
    --procs 1 2 4 \
    --solvers amg ddm \
    --output results_strong.csv \
    --benchmark-path ../cmake-build-release/benchmark/poisson_scaling/benchmark_poisson
```

#### Weak Scaling (Fixed Problem Per Processor)

```bash
cd benchmark
python3 run_scaling.py \
    --study weak \
    --base-mesh-size 32 \
    --order 2 \
    --procs 1 2 4 \
    --solvers amg ddm \
    --output results_weak.csv \
    --benchmark-path ../cmake-build-release/benchmark/poisson_scaling/benchmark_poisson
```

### 4. Generate Plots

First install required Python packages:

```bash
pip3 install matplotlib numpy
```

Then generate plots:

```bash
cd benchmark

# Strong scaling plots
python3 plot_results.py results_strong.csv \
    --study strong \
    --output strong_scaling.png

# Weak scaling plots
python3 plot_results.py results_weak.csv \
    --study weak \
    --output weak_scaling.png
```

## Benchmark Output Format

The benchmark application outputs parseable results:

```
=== Poisson Solver Benchmark ===
Solver type: ddm
Mesh size: 32x32
Polynomial order: 2
MPI processes: 4
Global DOFs: 4225
Using DDM solver

=== Results ===
DOFS: 4225
ASSEMBLE_TIME: 0.0160278
SOLVE_TIME: 0.0498474
TOTAL_TIME: 0.0658752
ITERATIONS: 33
L2_ERROR: 4.16612e-06
PROCS: 4
```

Key metrics:
- `DOFS`: Total degrees of freedom
- `ASSEMBLE_TIME`: Time to assemble linear system (seconds)
- `SOLVE_TIME`: Time to solve linear system (seconds)
- `TOTAL_TIME`: Total time (assemble + solve)
- `ITERATIONS`: Number of solver iterations
- `L2_ERROR`: L2 norm of solution error
- `PROCS`: Number of MPI processes

## Command-Line Options

### benchmark_poisson

```
Options:
  --solver TYPE      Solver type: 'amg' or 'ddm' (default: amg)
  --mesh-size N      Mesh elements per dimension (default: 16)
  --order N          Polynomial order (default: 2)
  --help             Print help message
```

### run_scaling.py

```
Options:
  --study TYPE           Scaling study: 'strong' or 'weak'
  --solvers SOLVERS      Solvers to benchmark: 'amg', 'ddm', or both (default: both)
  --mesh-size N          Mesh size for strong scaling (default: 64)
  --base-mesh-size N     Base mesh size for weak scaling (default: 32)
  --order N              Polynomial order (default: 2)
  --procs PROCS          Processor counts (default: 1 2 4)
  --output FILE          Output CSV file (required)
  --benchmark-path PATH  Path to benchmark executable
```

### plot_results.py

```
Options:
  input              Input CSV file from run_scaling.py
  --study TYPE       Scaling study: 'strong' or 'weak'
  --output FILE      Output plot file (PNG)
```

## Scaling Study Methodology

### Strong Scaling

Tests parallel efficiency with a **fixed problem size** and **varying number of processors**.

- **Ideal behavior**: Speedup = N (linear scaling)
- **Efficiency**: E = Speedup / N
- **Measures**: Communication overhead, load balancing

### Weak Scaling

Tests scalability with **fixed problem size per processor** and **increasing total problem size**.

- **Ideal behavior**: Time remains constant
- **Measures**: Algorithm scalability, memory bandwidth

## Example Results

Strong scaling with 64x64 mesh, order 2 (~16k DOFs):

| Solver | Procs | Time (s) | Speedup | Efficiency |
|--------|-------|----------|---------|------------|
| AMG    | 1     | 0.150    | 1.0     | 1.00       |
| AMG    | 2     | 0.082    | 1.83    | 0.91       |
| AMG    | 4     | 0.045    | 3.33    | 0.83       |
| DDM    | 1     | 0.180    | 1.0     | 1.00       |
| DDM    | 2     | 0.095    | 1.89    | 0.95       |
| DDM    | 4     | 0.055    | 3.27    | 0.82       |

## Troubleshooting

### "Benchmark executable not found"

Make sure to build the benchmark application first and provide the correct path:

```bash
cmake --build . --target benchmark_poisson
```

### MPI errors

Ensure MPI is properly installed and configured:

```bash
which mpirun
mpirun --version
```

### Python dependencies

Install required packages:

```bash
pip3 install matplotlib numpy
```

Or use a virtual environment:

```bash
python3 -m venv venv
source venv/bin/activate
pip install matplotlib numpy
```

## References

- Strong scaling: https://en.wikipedia.org/wiki/Scalability#Strong_scaling
- Weak scaling: https://en.wikipedia.org/wiki/Scalability#Weak_scaling
- Parallel efficiency: https://en.wikipedia.org/wiki/Parallel_efficiency
