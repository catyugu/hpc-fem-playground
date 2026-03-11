# 提示词2

* 现在耦合非线性求解的效率太低了，事实上你可以做如下操作：建立场与场之间的依赖图，然后如果A对B两个子图之间只有单向依赖，就可以先只对A中的场进行非线性迭代，收敛之后再用A中的场计算B，减少计算量。例如我们的例子中，热和电场有双向依赖，但是热膨胀位移场只对热场有依赖，所以可以先让热——电耦合收敛稳定，再利用温度场直接计算热膨胀位移。
* 移除一些没有什么意义的垃圾代码和无必要的抽象。
* 其他你认为值得改进的地方。应该使代码更简洁，更通用。
* 用以下代码验证结果，要求最大相对误差 < 1e-3
```bash
cd /home/catyugu/Projects/hpc-fem-playground &&
./cmake-build-release/example/busbar_pipeline cases/busbar/case.xml
python3 scripts/compare_comsol_results.py cases/busbar/result.txt results/busbar/mpfem_result.txt
```
