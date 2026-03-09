# 提示词2

* 引入MKL PARDISO稀疏求解器的支持，用于加速矩阵求解（作为新的线性求解器选项）
* 尽可能使用mfem原生的api而非手动的自由度数组访问
* 其他你认为值得改进的地方。
* 用以下代码验证结果，要求最大相对误差 <1%
```bash
cd /home/catyugu/Projects/hpc-fem-playground &&
./cmake-build-release/example/busbar_pipeline cases/busbar/case.xml
python3 scripts/compare_comsol_results.py cases/busbar/result.txt results/busbar/mpfem_result.txt
```