# Karhunen-Loève Transform
# Same as pca, this code is a demo on dataset 'ColorHistogram.asc'
import numpy as np

data = np.loadtxt('ColorHistogram.asc', dtype=np.float64)
data = data[:,1:]

means = np.mean(data, axis = 0)
data -= means

data = data.T
print(means.shape)

print(data.shape)

C = (data @ data.T) / data.shape[1]

eig, P = np.linalg.eig(C)
print("eig = ", eig)
print("P = ", P)

P1 = P.T[:16, :]
print("select first 16 eig: ", eig[:16])
#print(P.T @ C @ P)
#print(P1 @ C @ P1.T)

# compress and rebuild
data1 = P1.T @ (P1 @ data)

print("original data[0] vector: ", data[0])
print("rebuild  data[0] vector: ", data1[0])
