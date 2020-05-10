# Z-test for model compare

# input:
#          | model1 | model2 | ...
# dataset1 |  acc   |  acc   |
# dataset2 | ...
# ...

# output:
# model compare crosstable

import numpy as np
import math

def Z_test(acc, datan, threshold=1.96):
    n=acc.shape[0]
    result=np.zeros([n,n,3]) # win-loss-draw
    for i in range(n):
        for j in range(n):
            for k in range(acc.shape[1]):
                pa=acc[i,k]
                pb=acc[j,k]
                p=(pa+pb)/2
                Z=(pa-pb)/math.sqrt(2*p*(1-p)/datan[k])
                if Z>threshold:
                    result[i,j,0]+=1
                elif Z<-threshold:
                    result[i,j,1]+=1
                else:
                    result[i,j,2]+=1
            print("%d-%d-%d "%tuple(result[i,j,:]),end='')
        print()
    return result

if __name__=="__main__":
    acc=[
    [92.09,79.62,87.19],
    [85.51,76.81,84.78],
    [81.95,58.05,70.73],
    [95.14,95.99,96.42],
    [76.24,83.50,84.49],
    [85.80,77.54,85.07],
    [72.40,75.91,76.82],
    [70.90,74.70,74.40],
    [67.29,48.59,59.81],
    [80.00,84.07,83.70],
    [81.94,83.23,87.10],
    [85.33,78.80,72.61],
    [85.33,78.80,82.61],
    [94.67,95.33,96.00],
    [78.95,94.74,92.98],
    [73.34,73.16,73.56],
    [77.03,83.11,86.49],
    [74.35,76.04,76.95],
    [78.85,69.71,76.92],
    [83.72,70.04,98.33],
    [71.04,45.04,74.94],
    [94.38,96.63,98.88],
    [93.07,93.07,96.04],
    ]
    datan=[898,690,205,699,303,690,768,1000,214,270,155,368,351,150,57,3200,148,768,208,958,846,178,101]
    ret=Z_test(np.array(acc).T/100, np.array(datan))
    