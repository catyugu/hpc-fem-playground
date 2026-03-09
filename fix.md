# 提示词2

* 用以下代码验证结果，要求最大相对误差 <1%
```bash
cd /home/catyugu/Projects/hpc-fem-playground &&
./cmake-build-release/example/busbar_pipeline cases/busbar/case.xml
python3 scripts/compare_comsol_results.py cases/busbar/result.txt results/busbar/mpfem_result.txt
```
* 目前代码运行时有一定概率出错，输出如下，我怀疑可能是mfem内部多线程的问题？请修改我们的代码使其线程安全
```
(base) catyugu@LAPTOP-TLRCI986:~/Projects/hpc-fem-playground$ ./cmake-build-release/example/busbar_pipeline cases/busbar/case.xml
[INFO] XML parsing completed in 0.000s
[INFO] Material loading completed in 0.000s
[INFO] Problem building completed in 0.000s
[INFO] Reference result loading completed in 0.015s
Parsing cases/busbar/mesh.mphtxt ...
Parsed 7360 nodes and 39374 elements in 7 blocks (sdim=3)
Writing MFEM mesh to cases/busbar/mesh.mesh ...
Done.
[INFO] Mesh loading completed in 0.346s
[INFO] Using direct solver
[INFO] Using direct solver
[INFO] Using direct solver
[INFO] --- Coupling iteration 1 ---
[INFO] Electrostatics solve BC completed in 0.000s
[INFO] Electrostatics solve assemble completed in 0.025s
[INFO] Electrostatics solve solve completed in 0.160s
[INFO] Electrostatics solve breakdown: BC=0.000s, assemble=0.025s, solve=0.160s
[INFO] Electrostatics solve completed in 0.186s
[INFO] Heat transfer solve BC completed in 0.000s
[INFO] Heat transfer solve assemble completed in 0.144s
[ERROR] DirectSolver: RHS vector contains non-finite value at index 282 (value: inf)
(base) catyugu@LAPTOP-TLRCI986:~/Projects/hpc-fem-playground$ ./cmake-build-release/example/busbar_pipeline cases/busbar/case.xml
[INFO] XML parsing completed in 0.000s
[INFO] Material loading completed in 0.000s
[INFO] Problem building completed in 0.000s
[INFO] Reference result loading completed in 0.010s
Parsing cases/busbar/mesh.mphtxt ...
Parsed 7360 nodes and 39374 elements in 7 blocks (sdim=3)
Writing MFEM mesh to cases/busbar/mesh.mesh ...
Done.
[INFO] Mesh loading completed in 0.348s
[INFO] Using direct solver
[INFO] Using direct solver
[INFO] Using direct solver
[INFO] --- Coupling iteration 1 ---
[INFO] Electrostatics solve BC completed in 0.000s
[INFO] Electrostatics solve assemble completed in 0.024s
[INFO] Electrostatics solve solve completed in 0.160s
[INFO] Electrostatics solve breakdown: BC=0.000s, assemble=0.024s, solve=0.160s
[INFO] Electrostatics solve completed in 0.184s
[INFO] Heat transfer solve BC completed in 0.000s
[INFO] Heat transfer solve assemble completed in 0.147s
[INFO] Heat transfer solve solve completed in 0.149s
[INFO] Heat transfer solve breakdown: BC=0.000s, assemble=0.147s, solve=0.149s
[INFO] Heat transfer solve completed in 0.296s
[INFO] Solid mechanics solve BC completed in 0.000s
[INFO] Solid mechanics solve assemble completed in 0.219s
[INFO] Solid mechanics solve solve completed in 0.601s
[INFO] Solid mechanics solve breakdown: BC=0.000s, assemble=0.219s, solve=0.601s
[INFO] Solid mechanics solve completed in 0.820s
[INFO] --- Coupling iteration 2 ---
[INFO] Electrostatics solve BC completed in 0.000s
[ERROR] ConductivityCoefficient: resistivity is invalid (rho0=0.000000, alpha=0.003900, T=-67695839882836505982629483360470242094939500154348686229027919521657567269017378314282993111475937083392.000000, tref=298.000000, rho=-4541036939340672174664194518914899974650939072197581093769331274139395285401299310558840881152.000000)
```
* 其他你认为值得改进的地方。应该使代码更简洁，更通用。