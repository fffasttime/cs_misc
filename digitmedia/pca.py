from PIL import Image
import numpy as np
from matplotlib import pyplot as plt

BLOCK_SIZE = 4
PCA_DIM = 1

img = Image.open('test_B.bmp') # color image
img = np.array(img).astype(np.float32) / 256


print(img.shape, img.dtype)
row,col=img.shape[:2]

data=[]

for i in range(0,row,BLOCK_SIZE):
    for j in range(0,col,BLOCK_SIZE):
        data.append(img[i:i+BLOCK_SIZE,j:j+BLOCK_SIZE,:].reshape(-1))

data = np.array(data)
print(data.shape)

means = np.mean(data, axis = 0)
data -= means

data = data.T

#cov = data @ data.T
#diag, P = np.linalg.eig(cov)
#P=P.T

U, diag, V = np.linalg.svd(data, full_matrices=False)
P = U.T
print("eig value: ", diag)
P1 = P[:PCA_DIM,:]

output = P1 @ data

data1 = P1.T @ output
data1 = data1.T
data1 += means
print(data1.shape)
img1=np.zeros(img.shape)

for i in range(0,row,BLOCK_SIZE):
    for j in range(0,col,BLOCK_SIZE):
        img1[i:i+BLOCK_SIZE,j:j+BLOCK_SIZE,:]=\
            data1[(i//BLOCK_SIZE*(col//BLOCK_SIZE)+j//BLOCK_SIZE)]\
            .reshape((BLOCK_SIZE, BLOCK_SIZE, 3))

plt.subplot(121); plt.imshow(img)
plt.subplot(122); plt.imshow(img1)
plt.show()