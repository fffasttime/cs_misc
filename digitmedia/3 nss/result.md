

## 索引方法：iDistance索引或LSH索引

### LSH算法描述

LSH的主要思想是，高维空间的两点若距离很近，那么设计一种哈希函数对这两点进行哈希值计算，使得他们哈希值有很大的概率是一样的。同时若两点之间的距离较远，他们哈希值相同的概率会很小。距离敏感的哈希函数：哈希映射后的距离随所属概率单调。

最初的LSH使用L1代替L2范数，嵌入到Hamming空间。因为Hamming空间里有合适的哈希函数。

对于欧氏空间，可使用投影的方式。对点进行投影变换，如果两个点距离小，则投影后仍可能很小。因此可在具有相同哈希值的集合中寻找近邻点。如果单次查询准确率较低，可以使用多次哈希组合提高准确率。

代码实现中，construct()实现构建这样一个哈希表，首先它生成一系列随机向量用于投影。gethash()接受两个向量，计算两个向量的投影，得到哈希编码。

### 测试

在ColorHistogram数据集上执行前10近邻搜索1000次，并与暴力查找的准确结果进行比对计算准确率。

在使用长度为12的哈希编码，进行五次哈希组合的情况下，测得平均查找次数是2121.55，约为原始数据60840的三十分之一。

下表为对比几种参数下的测试结果。其中要求前10完全匹配准确率大概在60%左右。由于这是近似查询，完全匹配是一个很高的要求。可以看出相对于直接蛮力查询虽然牺牲了一定的准确率，在时间效率上还是很有优势的。


|round|length|top 10 completely match|top 10 recall|building time|find items time|
|--|--|--|--|--|--|
|r= 5|len= 12|acc=58.00%|recall=71.38%|421 ms|1778 ms|
|r= 5|len= 12|acc=58.00%|recall=71.38%| 742 ms| 2496 ms|
|r= 5|len= 12|acc=58.00%|recall=71.38%| 422 ms| 1676 ms|
|r= 5|len= 12|acc=58.00%|recall=71.38%| 413 ms| 1764 ms|
|r= 4|len= 10|acc=59.20%|recall=72.38%| 420 ms| 2837 ms|
|r= 6|len= 12|acc=63.70%|recall=76.60%| 510 ms| 2247 ms|
|r= 4|len= 11|acc=55.20%|recall=69.05%| 606 ms| 2340 ms|
|r= 5|len= 12|acc=58.00%|recall=71.38%| 684 ms| 2318 ms|