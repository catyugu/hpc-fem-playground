# 提示词2

* 现在的直接求解器对二阶而言效率太差了，请引入PARDISO求解器作为可选项（选择性编译）：
```text
(base) catyugu@LAPTOP-TLRCI986:~/Projects/hpc-fem-playground$ ./cmake-build-release/example/busbar_pipeline cases/busbar_order2/case.xml
[INFO] XML parsing completed in 0.000s
[INFO] Material loading completed in 0.000s
[INFO] Problem building completed in 0.000s
[INFO] Reference result loading completed in 0.010s
Parsing cases/busbar_order2/mesh.mphtxt ...
Parsed 7360 nodes and 39374 elements in 7 blocks (sdim=3)
Writing MFEM mesh to cases/busbar_order2/mesh.mesh ...
Done.
[INFO] Mesh loading completed in 0.299s
[INFO] Using direct solver
[INFO] Using direct solver
[INFO] Using direct solver
[INFO] === Phase 1: Tightly coupled iteration ===
[INFO] --- Coupling iteration 1 ---
[INFO] Electrostatics solve BC completed in 0.000s
[INFO] Electrostatics solve assemble completed in 0.148s
[INFO] Electrostatics solve solve completed in 1.586s
[INFO] Electrostatics solve breakdown: BC=0.000s, assemble=0.148s, solve=1.586s
[INFO] Electrostatics solve completed in 1.734s
[INFO] Heat transfer solve BC completed in 0.000s
[INFO] Heat transfer solve assemble completed in 0.423s
[INFO] Heat transfer solve solve completed in 1.388s
[INFO] Heat transfer solve breakdown: BC=0.000s, assemble=0.423s, solve=1.388s
[INFO] Heat transfer solve completed in 1.810s
[INFO] --- Coupling iteration 2 ---
[INFO] Electrostatics solve BC completed in 0.000s
[INFO] Electrostatics solve assemble completed in 0.168s
[INFO] Electrostatics solve solve completed in 1.471s
[INFO] Electrostatics solve breakdown: BC=0.000s, assemble=0.168s, solve=1.471s
[INFO] Electrostatics solve completed in 1.639s
[INFO] Heat transfer solve BC completed in 0.000s
[INFO] Heat transfer solve assemble completed in 0.465s
[INFO] Heat transfer solve solve completed in 1.571s
[INFO] Heat transfer solve breakdown: BC=0.000s, assemble=0.465s, solve=1.571s
[INFO] Heat transfer solve completed in 2.036s
[INFO] --- Coupling iteration 3 ---
[INFO] Electrostatics solve BC completed in 0.000s
[INFO] Electrostatics solve assemble completed in 0.180s
[INFO] Electrostatics solve solve completed in 1.379s
[INFO] Electrostatics solve breakdown: BC=0.000s, assemble=0.180s, solve=1.379s
[INFO] Electrostatics solve completed in 1.559s
[INFO] Heat transfer solve BC completed in 0.000s
[INFO] Heat transfer solve assemble completed in 0.433s
[INFO] Heat transfer solve solve completed in 1.435s
[INFO] Heat transfer solve breakdown: BC=0.000s, assemble=0.433s, solve=1.435s
[INFO] Heat transfer solve completed in 1.868s
[INFO] --- Coupling iteration 4 ---
[INFO] Electrostatics solve BC completed in 0.000s
[INFO] Electrostatics solve assemble completed in 0.169s
[INFO] Electrostatics solve solve completed in 1.397s
[INFO] Electrostatics solve breakdown: BC=0.000s, assemble=0.169s, solve=1.397s
[INFO] Electrostatics solve completed in 1.567s
[INFO] Heat transfer solve BC completed in 0.000s
[INFO] Heat transfer solve assemble completed in 0.445s
[INFO] Heat transfer solve solve completed in 1.483s
[INFO] Heat transfer solve breakdown: BC=0.000s, assemble=0.445s, solve=1.483s
[INFO] Heat transfer solve completed in 1.928s
[INFO] --- Coupling iteration 5 ---
[INFO] Electrostatics solve BC completed in 0.000s
[INFO] Electrostatics solve assemble completed in 0.178s
[INFO] Electrostatics solve solve completed in 1.476s
[INFO] Electrostatics solve breakdown: BC=0.000s, assemble=0.178s, solve=1.476s
[INFO] Electrostatics solve completed in 1.655s
[INFO] Heat transfer solve BC completed in 0.000s
[INFO] Heat transfer solve assemble completed in 0.433s
[INFO] Heat transfer solve solve completed in 1.425s
[INFO] Heat transfer solve breakdown: BC=0.000s, assemble=0.433s, solve=1.425s
[INFO] Heat transfer solve completed in 1.858s
[INFO] Tightly coupled fields converged after 5 iterations
[INFO] === Phase 2: Solving downstream fields ===
[INFO] Solid mechanics solve BC completed in 0.000s
[INFO] Solid mechanics solve assemble completed in 0.885s
[INFO] Solid mechanics solve solve completed in 7.514s
[INFO] Solid mechanics solve breakdown: BC=0.000s, assemble=0.885s, solve=7.514s
[INFO] Solid mechanics solve completed in 8.399s
[INFO] Coupled solve completed in 26.056s
[INFO] Result export completed in 0.021s
```
* 移除一些没有什么意义的垃圾代码和无必要的抽象。
* 其他你认为值得改进的地方。应该使代码更简洁，更通用。
* 用以下代码验证结果，要求最大相对误差 < 1e-3
```bash
cd /home/catyugu/Projects/hpc-fem-playground &&
./cmake-build-release/example/busbar_pipeline cases/busbar/case.xml
python3 scripts/compare_comsol_results.py cases/busbar/result.txt results/busbar/mpfem_result.txt
```
