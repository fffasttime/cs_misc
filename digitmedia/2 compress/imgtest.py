from PIL import Image
import numpy as np
from matplotlib import pyplot as plt

BLOCK_SIZE = 4

### read image ###

img = Image.open('test_B.bmp') # color image
img = np.array(img).astype(np.float32) / 256

print(img.shape, img.dtype)
row,col=img.shape[:2]

data=[]

for i in range(0,row,BLOCK_SIZE):
    for j in range(0,col,BLOCK_SIZE):
        data.append(img[i:i+BLOCK_SIZE,j:j+BLOCK_SIZE,:].reshape(-1))

data = np.array(data)

### read image ###

def rebuildimage(data1):
    img1=np.zeros(img.shape)

    for i in range(0,row,BLOCK_SIZE):
        for j in range(0,col,BLOCK_SIZE):
            img1[i:i+BLOCK_SIZE,j:j+BLOCK_SIZE,:]=\
                data1[(i//BLOCK_SIZE*(col//BLOCK_SIZE)+j//BLOCK_SIZE)]\
                .reshape((BLOCK_SIZE, BLOCK_SIZE, 3))

    return img1

import lbg

plt.subplot(231); plt.imshow(img); plt.title('input')

data1=lbg.lbg(data, 11)
img1=rebuildimage(data1)
plt.subplot(232); plt.imshow(img1); plt.title('lgb, 1/4')

data1=lbg.lbg(data, 10)
img1=rebuildimage(data1)
plt.subplot(233); plt.imshow(img1); plt.title('lgb, 1/8')

import pca

data1=pca.pca(data.copy(), 1)
img1=rebuildimage(data1)
plt.subplot(236); plt.imshow(img1); plt.title('klt, 1 eig, 1/48')

data1=pca.pca(data.copy(), 4)
img1=rebuildimage(data1)
plt.subplot(235); plt.imshow(img1); plt.title('klt, 4 eig, 1/12')

data1=pca.pca(data.copy(), 12)
img1=rebuildimage(data1)
plt.subplot(234); plt.imshow(img1); plt.title('klt, 12 eig, 1/4')

plt.show()