# Linde-Buzo-Gray Algorithm
# Vector Quantization

import numpy as np

def split(codevec, eps):
    ret=[]
    for v in codevec:
        ret.append(v + eps)
        ret.append(v - eps)
    return ret

def lbg(data, steps):
    print(data.shape)
    n=data.shape[0]
    m=data.shape[1]
    codevec=[np.mean(data, axis=0)]   
    coded=[0 for i in range(n)]
    for i in range(steps):
        codevec=split(codevec, np.zeros(m)+0.01)
        sumvec=[np.zeros(m) for x in range(len(codevec))]
        cntvec=[0 for x in range(len(codevec))]
        for j in range(n):
            dis1=np.linalg.norm(data[j]-codevec[coded[j]*2])
            dis2=np.linalg.norm(data[j]-codevec[coded[j]*2+1])
            if dis1<dis2:
                coded[j]*=2
            else:
                coded[j]=coded[j]*2+1
            sumvec[coded[j]]+=data[j]
            cntvec[coded[j]]+=1
        for j in range(len(codevec)):
            if cntvec[j]:
                codevec[j]=sumvec[j]/cntvec[j]
    print(data[1])
    print(codevec[coded[1]])
    ret=np.zeros((n,m))
    for i in range(n):
        ret[i,:]=codevec[coded[i]]
    return ret