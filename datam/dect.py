# a slow and poor decision tree

import numpy as np
import math

def infoEnpt(x):
    return -sum(map(lambda a:a*math.log2(a) if a else 0,x))
def infoEnpt2(x):
    if x==1 or x==0:
        return 0
    return -x*math.log2(x)-(1-x)*math.log2(1-x)

def gini(x):
    return 1-sum(map(lambda a:a*a,x))
def gini2(x):
    return 1-x**2-(1-x)**2

def minpb2(x):
    return min(x,1-x)

def pureness(x, f):
    n=sum(x)
    return f(map(lambda a:a/n,x))

# average pureness of children
def sum_pureness(x, f=infoEnpt):
    n=0
    s=0.0
    for x1 in x:
        n+=sum(x1)
        s+=sum(x1)*pureness(x1, f)
    return s/n

# average pureness of two children
# s1: positve, c1: first part, pos: positive in first part
def sum_pureness2(n, s1, c1, pos1, f=infoEnpt2):
    return (f(pos1/c1)*c1+f((s1-pos1)/(n-c1))*(n-c1))/n

# try spilt point for continous 
# try_value([3.0,1.0,2.0],[1,0,1])
def try_value(val, label, f=infoEnpt2):
    v1=sorted(zip(val, label))
    val,label=zip(*v1)
    n=len(val)
    s1=sum(label)
    pos=0
    for i in range(n-1):
        pos+=label[i]
        if i and val[i]==val[i-1]: 
            continue
        s=sum_pureness2(n,s1,i+1,pos)
        print(("<=%f : %f")%(val[i], s))

class MetaData:
    # attr: None for continious value
    def __init__(self, xname, attr):
        self.xname=xname
        self.attr=attr
        self.nx=len(xname)

class Node:
    def __init__(self):
        self.ch=[]
    def print(self, deep=0):
        for i in range(deep):
            print('  ',end='')
        print("%s : %d/%d"%(self.name,self.pos,self.n))
        for c in self.ch:
            c.print(deep+1)
    def fit(self, tree,datax,datay,freq,nodename,deep=0):
        self.name=nodename
        self.n=sum(freq)
        self.pos=sum(freq[datay==1])
        if deep==2: return
        if self.pos==0 or self.pos==self.n:
            return
        minsp=None ; mini=None
        curpure=sum_pureness([[self.pos,self.n-self.pos]], tree.mode)
        for i,x in enumerate(datax):
            if not tree.idlex[i]:
                continue
            if type(tree.metadata.attr[i]) is int:
                cnt=np.zeros((tree.metadata.attr[i],2))
                # print(datax[i])
                for j in range(len(freq)):
                    cnt[datax[i][j],datay[j]]+=freq[j]
                print(cnt)
                s=sum_pureness(cnt, tree.mode)
                print(nodename,':',tree.metadata.xname[i],s,curpure-s)
                if minsp==None or minsp>s:
                    minsp=s
                    mini=i
            else:
                pass
                # Todo : continious value
        if deep==0: mini=0
        if mini is not None:
            tree.idlex[mini]=False
            self.pid=mini
            for j in range(tree.metadata.attr[mini]):
                self.ch.append(Node())
                sel=datax[mini]==j
                newx=[x[sel] for x in datax]
                self.ch[-1].fit(tree, newx, datay[sel], freq[sel], "%s_%d"%(tree.metadata.xname[mini],j), deep+1)
                self.ch[-1].id=j

            tree.idlex[mini]=True

class Tree:
    def __init__(self, metadata, mode=infoEnpt):
        self.metadata=metadata
        self.mode=mode
    def fit(self, datax, datay, freq=None):
        self.root=Node()
        self.nitem=len(datax)
        self.idlex=np.ones(self.metadata.nx).astype(bool)
        if freq is None:
            freq=np.ones(self.nitem)
        self.root.fit(self,datax,datay, freq,'root')
    def print(self):
        self.root.print()
    def predict(self, datax):
        node=self.root
        while len(node.ch):
            node=node.ch[datax[node.pid]]
        return node.pos>node.n-node.pos
            

if __name__=="__main__":
    datax=np.array([
        [0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1],
        [0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1],
        [0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1],
    ])
    datay=np.array([1,0]*8)
    freq=np.array([5,40,0,15,10,5,45,0,10,5,25,0,5,20,0,15])
    freq=np.array([5,0,0,20,20,0,0,5,0,0,25,0,0,0,0,25])

    tree=Tree(MetaData("X Y Z".split(),[2,2,2]),min)
    tree.fit(datax, datay, freq)
    tree.print()

    acc=0
    for i in range(len(datay)):
        ret=tree.predict([x[i] for x in datax])
        if ret==datay[i]:
            acc+=freq[i]
    print(acc/sum(freq))