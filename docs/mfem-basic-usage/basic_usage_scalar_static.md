# Basic usage — scalar static (steady) field with MFEM

This short guide shows the minimal pattern to solve steady-state scalar diffusion/heat problems with MFEM. It is distilled from the two examples in this repo:

- `example/ex0_cylinder_dirichlet.cpp` — purely Dirichlet example.
- `example/ex1_cube_mixed_bc.cpp` — Dirichlet + Robin + Neumann mixed example.

Keep this as a quick recipe you can follow when you need a simple scalar Poisson/heat solve.

## Contract (inputs / outputs)

- Inputs: a mesh file (Gmsh or MFEM mesh) and boundary attribute markers.
- Output: a GridFunction containing the scalar field (saved with `ParaViewDataCollection`).

## Minimal steps (2–10 lines each)

1) Read the mesh and create the H1 finite element space

   - Mesh: `Mesh mesh(mesh_file);`
   - Space: `H1_FECollection fec(order, dim); FiniteElementSpace fes(&mesh, &fec);`

2) Define the BilinearForm (left-hand side)

   - Diffusion term: `BilinearForm a(&fes); a.AddDomainIntegrator(new DiffusionIntegrator(k_coef));`
   - If you have Robin (convective) BCs of form `-k dT/dn = h (T - T_inf)`, add a boundary mass term
     `a.AddBoundaryIntegrator(new MassIntegrator(h_coef), robin_bdr);`

3) Define Dirichlet (essential) BCs

   - Choose attributes that are essential (e.g., cold boundary). Build `Array<int> ess_bdr` and call
     `fes.GetEssentialTrueDofs(ess_bdr, ess_tdof_list);`
   - Project known Dirichlet values into the `GridFunction x` using
     `x.ProjectBdrCoefficient(dirichlet_coef, ess_bdr);`

4) Build the LinearForm (right-hand side)

   - For Robin RHS contribution: `lf.AddBoundaryIntegrator(new BoundaryLFIntegrator(h*T_inf), robin_bdr);`
   - For Neumann (prescribed flux q): `lf.AddBoundaryIntegrator(new BoundaryLFIntegrator(q_coef), neumann_bdr);`

5) Assemble and form the linear system

   - `a.Assemble(); lf.Assemble();`
   - `a.FormLinearSystem(ess_tdof_list, x, lf, A, X, B);` — this enforces Dirichlet BCs in the linear system.

6) Solve

   - Use any MFEM solver; examples use `CGSolver` with a simple preconditioner: `GSSmoother`.

7) Recover the full FEM solution

   - `a.RecoverFEMSolution(X, lf, x);`

8) Save / postprocess

   - Use `ParaViewDataCollection` to write `x` for visualization.

## Notes and common pitfalls

- Sign convention: `BoundaryLFIntegrator(q)` adds \int_Gamma q v to the RHS. In the examples a positive `q` means heat entering the domain. If you prefer the opposite sign, pass `-q`.
- Robin BC weak form used here adds `\int_Gamma h u v` to the left-hand side and `\int_Gamma h T_inf v` to the RHS.
- Always assemble the bilinear form (`a.Assemble()`) and linear form (`lf.Assemble()`) after adding all domain and boundary integrators.
- Use `FormLinearSystem` with the assembled `LinearForm` to correctly include natural boundary contributions and enforce essential DOFs.

## Small example mapping

- Pure Dirichlet (see `ex0_cylinder_dirichlet.cpp`): build `a` with diffusion, set `ess_bdr`, project Dirichlet values into `x`, then `a.FormLinearSystem(ess_tdof_list, x, b, A, X, B)` where `b` is a zero vector or domain load.
- Mixed (see `ex1_cube_mixed_bc.cpp`): add `MassIntegrator(h)` on the Robin boundary, add `BoundaryLFIntegrator(h*T_inf)` for the Robin RHS and `BoundaryLFIntegrator(q)` for Neumann flux, then `a.FormLinearSystem(ess_tdof_list, x, lf, A, X, B)`.

## Quick commands (build & run)

Use the repository build from this workspace (example target):

```bash
cmake --build cmake-build-debug --target example -j 4
./cmake-build-debug/example/ex1_cube_mixed_bc [mesh_file] [output_dir]
```

## Next steps / extensions

- Make boundary parameters (`h`, `T_inf`, `q`) configurable via CLI or input file.
- Add time-dependence by using a time integrator and mass matrix (Transient heat).
- Replace `CGSolver` + `GSSmoother` with Hypre or AMG for larger problems.

This file aims to be the minimal reference for the simplest scalar static problems in MFEM. If you'd like, I can add a tiny template C++ file that implements this recipe with options parsing and unit tests.
