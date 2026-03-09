# 提示词2

* mpi并行存在问题，数据交换不稳定，即使是没有开启np选项，有时候也会出现不合理的数据。
* 问题可能出在数据的共享与通信上。你可以参考mfem的并行示例和应用来寻找问题所在。你可能应该用rank 0统一管理进程间的数据。
* 尽可能使用mfem原生的api而非手动的自由度数组访问
* 其他你认为值得改进的地方。
* 用以下代码验证结果，要求最大相对误差 <1%
```bash
cd /home/catyugu/Projects/hpc-fem-playground &&
mpirun -np 1 ./cmake-build-release/example/busbar_pipeline cases/busbar/case.xml
python3 scripts/compare_comsol_results.py cases/busbar/result.txt results/busbar/mpfem_result.txt
```
* 验收标准：MPI（1个、2个、4个并行）运行至少连续三次不崩溃，结果正确