import numpy as np
from matplotlib import pyplot as plt

data=np.loadtxt('out.txt')

plt.hist(data, edgecolor='black',alpha=0.7)

plt.xlabel('result')
plt.ylabel('count')
plt.show()