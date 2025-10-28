#!/usr/bin/env python3
"""
Scaling study automation script for Poisson solver benchmarks.

This script runs strong and weak scaling studies:
- Strong scaling: Fixed problem size, varying number of processors
- Weak scaling: Fixed problem per processor, increasing total problem size

Usage:
    python3 run_scaling.py --study strong --output results_strong.csv
    python3 run_scaling.py --study weak --output results_weak.csv
"""

import subprocess
import sys
import csv
import argparse
import re
import os

# Constants
DEFAULT_STRONG_MESH_SIZE = 64
DEFAULT_STRONG_ORDER = 2
DEFAULT_WEAK_BASE_SIZE = 32
DEFAULT_WEAK_ORDER = 2
DEFAULT_PROC_COUNTS = [1, 2, 4]
BENCHMARK_EXECUTABLE = "./benchmark_poisson"


def parse_benchmark_output(output_text):
    """
    Parse benchmark output and extract key metrics.
    
    Expected format:
        DOFS: 4225
        ASSEMBLE_TIME: 0.0160278
        SOLVE_TIME: 0.0498474
        TOTAL_TIME: 0.0658752
        ITERATIONS: 33
        L2_ERROR: 4.16612e-06
        PROCS: 4
    
    Returns:
        dict with keys: dofs, assemble_time, solve_time, total_time, iterations, l2_error, procs
        Returns None if parsing fails
    """
    result = {}
    
    patterns = {
        'dofs': r'DOFS:\s+(\d+)',
        'assemble_time': r'ASSEMBLE_TIME:\s+([\d.e+-]+)',
        'solve_time': r'SOLVE_TIME:\s+([\d.e+-]+)',
        'total_time': r'TOTAL_TIME:\s+([\d.e+-]+)',
        'iterations': r'ITERATIONS:\s+(\d+)',
        'l2_error': r'L2_ERROR:\s+([\d.e+-]+)',
        'procs': r'PROCS:\s+(\d+)'
    }
    
    for key, pattern in patterns.items():
        match = re.search(pattern, output_text)
        if match:
            if key in ['dofs', 'iterations', 'procs']:
                result[key] = int(match.group(1))
            else:
                result[key] = float(match.group(1))
        else:
            print(f"Warning: Could not parse {key} from output", file=sys.stderr)
            return None
    
    return result


def run_benchmark(solver, mesh_size, order, num_procs, benchmark_path):
    """
    Run a single benchmark configuration.
    
    Returns:
        dict with benchmark results, or None if run failed
    """
    cmd = [benchmark_path, 
           "--solver", solver, 
           "--mesh-size", str(mesh_size), 
           "--order", str(order)]
    
    if num_procs > 1:
        cmd = ["mpirun", "-np", str(num_procs)] + cmd
    
    print(f"Running: {' '.join(cmd)}", file=sys.stderr)
    
    try:
        result = subprocess.run(cmd, 
                              capture_output=True, 
                              text=True, 
                              timeout=300,
                              check=True)
        output = result.stdout
        
        metrics = parse_benchmark_output(output)
        if metrics is None:
            print(f"Failed to parse output for configuration: "
                  f"solver={solver}, mesh_size={mesh_size}, "
                  f"order={order}, procs={num_procs}", file=sys.stderr)
            return None
        
        # Add configuration info
        metrics['solver'] = solver
        metrics['mesh_size'] = mesh_size
        metrics['order'] = order
        
        return metrics
        
    except subprocess.TimeoutExpired:
        print(f"Timeout for configuration: solver={solver}, "
              f"mesh_size={mesh_size}, order={order}, procs={num_procs}", 
              file=sys.stderr)
        return None
    except subprocess.CalledProcessError as e:
        print(f"Benchmark failed with exit code {e.returncode}", file=sys.stderr)
        print(f"stderr: {e.stderr}", file=sys.stderr)
        return None


def strong_scaling_study(solvers, mesh_size, order, proc_counts, output_file, benchmark_path):
    """
    Run strong scaling study: fixed problem size, varying processors.
    
    For each solver, run with different processor counts.
    """
    print(f"\n=== Strong Scaling Study ===", file=sys.stderr)
    print(f"Problem size: {mesh_size}x{mesh_size}, order {order}", file=sys.stderr)
    print(f"Processor counts: {proc_counts}", file=sys.stderr)
    print(f"Solvers: {solvers}", file=sys.stderr)
    
    results = []
    
    for solver in solvers:
        for num_procs in proc_counts:
            metrics = run_benchmark(solver, mesh_size, order, num_procs, benchmark_path)
            if metrics:
                results.append(metrics)
            else:
                print(f"Warning: Skipping failed run: "
                      f"solver={solver}, procs={num_procs}", file=sys.stderr)
    
    # Write to CSV
    if results:
        fieldnames = ['solver', 'mesh_size', 'order', 'procs', 'dofs', 
                     'assemble_time', 'solve_time', 'total_time', 
                     'iterations', 'l2_error']
        
        with open(output_file, 'w', newline='') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(results)
        
        print(f"\nResults written to {output_file}", file=sys.stderr)
        print(f"Total runs: {len(results)}", file=sys.stderr)
    else:
        print("ERROR: No successful runs", file=sys.stderr)
        sys.exit(1)


def weak_scaling_study(solvers, base_mesh_size, order, proc_counts, output_file, benchmark_path):
    """
    Run weak scaling study: fixed problem per processor.
    
    For each processor count, scale mesh size to maintain constant DOFs per processor.
    Approximate scaling: mesh_size proportional to sqrt(num_procs)
    """
    print(f"\n=== Weak Scaling Study ===", file=sys.stderr)
    print(f"Base mesh size (1 proc): {base_mesh_size}x{base_mesh_size}, order {order}", 
          file=sys.stderr)
    print(f"Processor counts: {proc_counts}", file=sys.stderr)
    print(f"Solvers: {solvers}", file=sys.stderr)
    
    results = []
    
    for solver in solvers:
        for num_procs in proc_counts:
            # Scale mesh size with sqrt(num_procs) to maintain DOFs/proc
            scaled_mesh_size = int(base_mesh_size * (num_procs ** 0.5))
            
            print(f"\nProcs: {num_procs}, scaled mesh: {scaled_mesh_size}x{scaled_mesh_size}", 
                  file=sys.stderr)
            
            metrics = run_benchmark(solver, scaled_mesh_size, order, num_procs, benchmark_path)
            if metrics:
                # Add DOFs per processor
                metrics['dofs_per_proc'] = metrics['dofs'] / num_procs
                results.append(metrics)
            else:
                print(f"Warning: Skipping failed run: "
                      f"solver={solver}, procs={num_procs}", file=sys.stderr)
    
    # Write to CSV
    if results:
        fieldnames = ['solver', 'mesh_size', 'order', 'procs', 'dofs', 'dofs_per_proc',
                     'assemble_time', 'solve_time', 'total_time', 
                     'iterations', 'l2_error']
        
        with open(output_file, 'w', newline='') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(results)
        
        print(f"\nResults written to {output_file}", file=sys.stderr)
        print(f"Total runs: {len(results)}", file=sys.stderr)
    else:
        print("ERROR: No successful runs", file=sys.stderr)
        sys.exit(1)


def main():
    parser = argparse.ArgumentParser(
        description='Run scaling studies for Poisson solver benchmarks')
    
    parser.add_argument('--study', 
                       choices=['strong', 'weak'],
                       required=True,
                       help='Type of scaling study')
    
    parser.add_argument('--solvers',
                       nargs='+',
                       default=['amg', 'ddm'],
                       choices=['amg', 'ddm'],
                       help='Solvers to benchmark (default: both)')
    
    parser.add_argument('--mesh-size',
                       type=int,
                       help=f'Mesh size for strong scaling (default: {DEFAULT_STRONG_MESH_SIZE})')
    
    parser.add_argument('--base-mesh-size',
                       type=int,
                       help=f'Base mesh size for weak scaling (default: {DEFAULT_WEAK_BASE_SIZE})')
    
    parser.add_argument('--order',
                       type=int,
                       default=DEFAULT_STRONG_ORDER,
                       help=f'Polynomial order (default: {DEFAULT_STRONG_ORDER})')
    
    parser.add_argument('--procs',
                       nargs='+',
                       type=int,
                       default=DEFAULT_PROC_COUNTS,
                       help=f'Processor counts (default: {DEFAULT_PROC_COUNTS})')
    
    parser.add_argument('--output',
                       required=True,
                       help='Output CSV file')
    
    parser.add_argument('--benchmark-path',
                       default=BENCHMARK_EXECUTABLE,
                       help=f'Path to benchmark executable (default: {BENCHMARK_EXECUTABLE})')
    
    args = parser.parse_args()
    
    # Check if executable exists
    if not os.path.isfile(args.benchmark_path):
        print(f"ERROR: Benchmark executable not found: {args.benchmark_path}", 
              file=sys.stderr)
        print(f"Please build the benchmark or specify --benchmark-path", 
              file=sys.stderr)
        sys.exit(1)
    
    # Run requested study
    if args.study == 'strong':
        mesh_size = args.mesh_size if args.mesh_size else DEFAULT_STRONG_MESH_SIZE
        strong_scaling_study(args.solvers, mesh_size, args.order, 
                            args.procs, args.output, args.benchmark_path)
    elif args.study == 'weak':
        base_mesh_size = args.base_mesh_size if args.base_mesh_size else DEFAULT_WEAK_BASE_SIZE
        weak_scaling_study(args.solvers, base_mesh_size, args.order, 
                          args.procs, args.output, args.benchmark_path)


if __name__ == '__main__':
    main()
