# 提示词2

* 现在我们的代码不够优雅，含有参数非常多的超级函数。
* 性能较低且瓶颈未知，可能需要更多带时间戳的日志信息或者外部工具来确定瓶颈，并进行性能改善。
* 物理场的离散多项式阶数应该取为2阶以与COMSOL中配置相符。
* 其他你认为值得改进的地方。
* 用以下代码验证结果，要求最大相对误差 <1%
```bash
cd /home/catyugu/Projects/hpc-fem-playground &&
./cmake-build-release/example/busbar_pipeline cases/busbar/case.xml 
python3 scripts/compare_comsol_results.py cases/busbar/result.txt results/busbar/mpfem_result.txt
```