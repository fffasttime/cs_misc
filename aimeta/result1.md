Rosenbrock function  
$$\max f(x_1,x_2)=100(x_2-x_1^2)^2+(1-x_1)^2$$
$$s.t. -2.048\leq x_1,x_2 \leq 2.048$$  
$$x_1^*=x_2^*=-2.048, \quad v^*\approx 3905$$

# Compare table
|Task: argmax|parameters|average speed|accuracy(~3905)|
|---|---|---|---|
|GA(遗传算法)|ITER=1000 N=10 E=2 PC=0.6 PM=0.1|11.24ms|50.60% (253/500)|
||ITER=500 N=50 E=2 PC=0.6 PM=0.1|26.32ms|58.40% (292/500)|
||ITER=500 N=100 E=2 PC=0.6 PM=0.1|44.58ms|61.80% (309/500)|
|GA_noselection|ITER=1000 N=10 PC=0.6 PM=0.1|5.61ms|52.40% (262/500)|
||ITER=500 N=50 PC=0.6 PM=0.1|28.58ms|56.00% (280/500)|
|SA(模拟退火算法)|delta=0.99|**0.26ms**|49.00% (245/500)|
||delta=0.995|0.39ms|50.20% (251/500)|
||delta=0.99 mc_round=5|2.21ms|98.00% (490/500)|
||delta=0.99 mc_round=10|4.38ms|**100.00% (500/500)**|

## 蒙特卡洛方法  
适当多次运行取最优解的方式是蒙特卡洛方法。
设单步运行时间为$t$，准确度为$p$，要使误差小于$\varepsilon$  
则期望运行时间为
$$T=t\times log_{1-p}{\varepsilon}$$
这说明在本问题中，适当使用蒙特卡洛在运行效率上是十分划算的。
