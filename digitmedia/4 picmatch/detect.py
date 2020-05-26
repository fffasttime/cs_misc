# for first image in test data,
# find most similar image in train data using SIFT
import cv2
import numpy as np
import os
import shutil
import mysift
from PIL import Image

pictype=1
magic=[3,6,25,0,22]
source_file='cifar-10/test/%d/%d_%d.jpg'%(pictype,pictype,magic[pictype])
shutil.copyfile(source_file, 'result/'+source_file.split('/')[-1])

#read image
img1 = cv2.imread(source_file, cv2.IMREAD_COLOR)
img1 = cv2.cvtColor(img1,cv2.COLOR_BGR2GRAY)

if 0:
    sift=cv2.xfeatures2d.SIFT_create(6)
else:
    sift=mysift

kp1, des1 = sift.detectAndCompute(img1, None)

dis_min=1e9

train_folder='cifar-10/train/%d/'%pictype
cnt=0

for path in os.listdir(train_folder):
    if path.split('.')[-1]!='jpg' : continue
    #path = "0_19503.jpg"
    img2 = cv2.imread(train_folder+path, cv2.IMREAD_COLOR)
    img2 = cv2.cvtColor(img2,cv2.COLOR_BGR2GRAY)

    # sift search 2
    kp2, des2 = sift.detectAndCompute(img2, None)
    if des2 is None : continue
    
    # brute force matching
    bf = cv2.BFMatcher_create()
    matches = bf.match(des1, des2)

    dis_sum = 0

    for matche in matches:
        dis_sum += matche.distance

    cnt+=1
    if cnt%50==0: print(path, dis_sum)
    if dis_sum<dis_min:
        dis_min=dis_sum
        pic_similar=path

print(pic_similar, dis_min)

shutil.copyfile(train_folder+pic_similar, 'result/'+pic_similar)