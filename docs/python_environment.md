# Python Environment Configuration

## Overview

This project uses a dedicated conda environment `hpc-fem-playground` for all Python scripts and tools.

## Setup

### Create the environment

```bash
conda create -n hpc-fem-playground python=3.10
conda activate hpc-fem-playground
```

### Install required packages

```bash
pip install matplotlib numpy
```

## Usage

### Always activate before running Python scripts

```bash
conda activate hpc-fem-playground
```

### Running benchmark scripts

```bash
# Activate environment first
conda activate hpc-fem-playground

# Run scaling study
cd benchmark
python3 run_scaling.py --study strong --output results.csv --benchmark-path ../cmake-build-release/benchmark/poisson_scaling/benchmark_poisson

# Generate plots
python3 plot_results.py results.csv --study strong --output plot.png
```

### Deactivate when done

```bash
conda deactivate
```

## VSCode Integration

Add to `.vscode/settings.json`:

```json
{
    "python.defaultInterpreterPath": "~/miniconda3/envs/hpc-fem-playground/bin/python",
    "python.terminal.activateEnvironment": true
}
```

## Verify environment

```bash
conda activate hpc-fem-playground
python -c "import matplotlib, numpy; print('OK')"
```

Expected output: `OK`

## Troubleshooting

### "ModuleNotFoundError: No module named 'matplotlib'"

Make sure the conda environment is activated:
```bash
conda activate hpc-fem-playground
pip install matplotlib numpy
```

### Wrong Python version

Check which Python is being used:
```bash
which python
# Should show: ~/miniconda3/envs/hpc-fem-playground/bin/python
```

If not, activate the environment:
```bash
conda activate hpc-fem-playground
```
