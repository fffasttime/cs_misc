# Karhunen-Loève Transform

import numpy as np

data = np.loadtxt('ColorHistogram.asc', dtype=np.float64)
data = data[:,1:]

data = data.T

data = data - np.mean(data, axis = 0)

print(data.shape)

C = (data @ data.T) / data.shape[0]

eig, P = np.linalg.eig(C)
print("eig = ", eig)
print("P = ", P)

print("select first 10 eig: ", eig[:10])

P1 = P[:,:10]

print(P.T @ C @ P)
print(P1.T @ C @ P1)

print("original data[0] vector: ", data[0])
print("rebuild  data[0] vector: ", (P1.T@data)[0])

print((P1.T@data).shape)