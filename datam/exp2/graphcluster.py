import numpy as np

class Graph:
    def __init__(self):
        self.g=np.zeros([50,50],dtype=int)
        self.m=0
        self.deg=np.zeros([50],dtype=int)
        self.nodec=0

    def load(self, filename):
        with open(filename,'r') as f:
            print(f.readline())
            f.readline() #graph
            f.readline() #[
            while 1:
                s=f.readline().strip()
                if s==']':
                    break
                f.readline() #[
                source=f.readline().strip()
                if s=='node':
                    self.nodec+=1
                else:
                    target=f.readline().strip()
                    source=int(source.split(' ')[1])
                    target=int(target.split(' ')[1])
                    self.g[source,target]=1
                    self.g[target,source]=1
                    self.m+=1
                    self.deg[source]+=1
                    self.deg[target]+=1
                f.readline() #]

    def modulate(self, groups):
        sum=0.0
        for s in groups:
            for x in s:
                for y in s:
                    if x!=y:
                        sum+=self.g[x,y]-self.deg[x]*self.deg[y]/2/self.m
        return sum/self.m

    def modinc(self, g1, g2):
        sum=0.0
        for x in g1:
            for y in g2:
                sum+=self.g[x,y]-self.deg[x]*self.deg[y]/2/self.m
        return sum/2/self.m

def run():
    g=Graph()
    g.load('karate.gml')
    groups=[[i+1] for i in range(g.nodec)]
    for step in range(g.nodec-1):
        res=(0, 0, 0)
        for i in range(len(groups)):
            for j in range(i+1,len(groups)):
                inc=g.modinc(groups[i],groups[j])
                if inc>res[0]:
                    res=(inc,i,j)
        print(res)
        if res[0]<=0:
            break
        groups[res[1]]+=groups[res[2]]
        groups.pop(res[2])
        print(g.modulate(groups))
    print(groups)

if __name__=="__main__":
    run()
