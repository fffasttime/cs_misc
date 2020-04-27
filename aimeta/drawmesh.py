import numpy as np
import math
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import pyplot as plt

#XMIN=-3.0; YMIN=4.1
XMIN=10.1; YMIN=5.4
XMAX=12.1; YMAX=5.8
plt.xlim((XMIN-0.2,XMAX+0.2))
plt.ylim((YMIN,YMAX))
X=np.linspace(XMIN, XMAX,100)
Y=np.linspace(YMIN, YMAX,100)
X,Y=np.meshgrid(X,Y)
f=lambda x,y:21.5 + x*np.sin(4*math.pi*x) + y*np.sin(20*math.pi*y)
Z=f(X,Y)

plt.contourf(X,Y,Z)
plt.colorbar()

s=np.loadtxt('out.txt').reshape((-1,2))
#clustering
clust=[]
clustc=[]
for i in range(s.shape[0]):
    for j in range(len(clust)):
        if np.linalg.norm(s[i]-clust[j])<0.05:
            clustc[j]+=1
            break
    else:
        clust.append(s[i])
        clustc.append(1)
clust=np.array(clust)
plt.scatter(clust[:,0],clust[:,1], np.array(clustc), marker='x');
for i in range(clust.shape[0]):
    plt.annotate('%.2f,%.2f\n%.4f\n\n%d'%(clust[i,0],clust[i,1], f(clust[i,0],clust[i,1]) ,clustc[i]), \
        xy=(clust[i,0],clust[i,1]),xytext = (0, 28), textcoords = 'offset points',ha = 'center', va = 'top')

plt.show()
