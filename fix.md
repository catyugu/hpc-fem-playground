# 提示词2

* 我们未来要为每个物理场支持不同类型的边界条件，每种类型可能对应不同的边界编号，对应不同的系数。现在物理场接口的设计显然不够通用。同时也需要改变xml解析逻辑以支持更灵活、可扩展的边界条件施加。为了测试，尝试把case.xml中的某些边界条件拆分成多组，验证结果是不是一致。
* 移除一些没有什么意义的垃圾代码和过度抽象。
* 其他你认为值得改进的地方。应该使代码更简洁，更通用。
* 用以下代码验证结果，要求最大相对误差 < 1e-3
```bash
cd /home/catyugu/Projects/hpc-fem-playground &&
./cmake-build-release/example/busbar_pipeline cases/busbar/case.xml
python3 scripts/compare_comsol_results.py cases/busbar/result.txt results/busbar/mpfem_result.txt
```
