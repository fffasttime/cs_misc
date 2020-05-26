import cv2
import numpy as np
import math

def detectAndCompute(img, _):
    keypoint=cv2.goodFeaturesToTrack(img, 6, 0.01, 5).astype(int)
    #print(keypoint)
    img=img[:]
    img=cv2.GaussianBlur(img, (5,5),1,1)

    # TODO: Gauss pyramid

    # grad filter
    kernel = np.array([
        [[-1,0,1],[-1,0,1],[-1,0,1]],
        [[-1,-1,-1,],[0,0,0],[1,1,1]]
    ])
    gx=cv2.filter2D(img, -1, kernel[1])
    gy=cv2.filter2D(img, -1, kernel[0])
    grad=(gx**2+gy**2)**0.5
    angle=np.arctan2(gy, gx)

    # TODO: rotate

    desc=np.zeros([6,3,3,8],dtype=np.float32)
    angle_dir8=((1-angle/math.pi)*4).astype(int)

    for c, kp in enumerate(keypoint):
        x=kp[0,0]-4
        y=kp[0,1]-4
        for i in range(9):
            for j in range(9):
                xx=x+i
                yy=y+j
                if xx>=0 and yy>=0 and xx<32 and yy <32:
                    desc[c,i//3,j//3,angle_dir8[xx,yy]]+=1

    return keypoint, desc.reshape([6,72])