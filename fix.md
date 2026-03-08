# 提示词2

* 现在我们的代码不够优雅，尤其是mfem_coupled_solver.cpp，将各种求解都耦合在一起，而且含有参数非常多的超级函数。
* 我们需要分离各个物理场的定义、边界施加和源项设置。
* 需要改变下游代码和材料库互动的方式，使代码更加简洁和灵活。
* 物理场的离散多项式阶数不应该硬编码，而应该分别编码在配置文件的不同物理场中。
* 线性求解需要提供多种不同策略以引入更高效的求解器等。建议使用HYPRE求解器以进一步加速。
* 非线性求解应该使用牛顿-拉普松残差迭代方案。
* 其他你认为值得改进的地方。
* 用以下代码验证结果，要求最大相对误差 <2%
```bash
cd /home/catyugu/Projects/hpc-fem-playground &&
./cmake-build-release/example/busbar_pipeline
   cases/busbar/case.xml 
   python3 scripts/compare_comsol_results.py 
   cases/busbar/result.txt 
   results/busbar/mpfem_result.txt
```