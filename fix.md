# 提示词2

* 用以下代码验证结果，要求最大相对误差 <1%
```bash
cd /home/catyugu/Projects/hpc-fem-playground &&
./cmake-build-release/example/busbar_pipeline cases/busbar/case.xml
python3 scripts/compare_comsol_results.py cases/busbar/result.txt results/busbar/mpfem_result.txt
```
* 目前代码运行时概率出错，输出如下：
```
(base) catyugu@LAPTOP-TLRCI986:~/Projects/hpc-fem-playground$ ./cmake-build-release/example/busbar_pipeline cases/busbar/case.xml
[INFO] XML parsing completed in 0.000s
[INFO] Material loading completed in 0.000s
[INFO] Problem building completed in 0.000s
[INFO] Reference result loading completed in 0.010s
Parsing cases/busbar/mesh.mphtxt ...
Parsed 7360 nodes and 39374 elements in 7 blocks (sdim=3)
Writing MFEM mesh to cases/busbar/mesh.mesh ...
Done.
[INFO] Mesh loading completed in 0.326s
[INFO] --- Coupling iteration 1 ---
[INFO] Electrostatics solve BC completed in 0.000s
[INFO] Electrostatics solve assemble completed in 0.023s
[INFO] Electrostatics solve solve completed in 0.060s
[INFO] Electrostatics solve breakdown: BC=0.000s, assemble=0.023s, solve=0.060s
[INFO] Electrostatics solve completed in 0.083s
[INFO] Heat transfer solve BC completed in 0.000s
[INFO] Heat transfer solve assemble completed in 0.077s


Verification failed: (IsFinite(nom)) is false:
 --> nom = -inf
 ... in function: virtual void mfem::CGSolver::Mult(const mfem::Vector&, mfem::Vector&) const
 ... in file: /home/catyugu/Projects/hpc-fem-playground/share/mfem/linalg/solvers.cpp:779

Aborted (core dumped)
```
* 其他你认为值得改进的地方。应该使代码更简洁，更通用。