#include <iostream>
#include <cstdio>
#include <cmath>
#include <cstring>
#include <ctime>
#include <algorithm>
using namespace std;

const int N=100010;
const double inf=1e10;
int n,m,K,rt;
double ans[2]; int ansp[2];
#define inc(i,n) for (int i=0;i<n;i++)

//s[]:tree node  p[2]:2-d coord of leaf node  x[2]:min(LB) coord of a subspace  y[2]:max(RT) coord
struct Node{
	double p[2],x[2],y[2]; int id;
	bool operator<(const Node &v)const{
		return p[K]<v.p[K];
	}
}a[N],s[N],q;
int ch[N][2];
#define lc ch[u][0]
#define rc ch[u][1]
void upd(int u){
	inc(i,2){
		if (lc) s[u].x[i]=min(s[u].x[i],s[lc].x[i]),
				s[u].y[i]=max(s[u].y[i],s[lc].y[i]);
		if (rc) s[u].x[i]=min(s[u].x[i],s[rc].x[i]),
				s[u].y[i]=max(s[u].y[i],s[rc].y[i]);
	}
}
void add(int u, Node &t){
	inc(i,2) s[u].x[i]=s[u].y[i]=s[u].p[i]=t.p[i];
	s[u].id=t.id;
}
double disL1Min(int u, Node &t){ //min L1 dis to a Rect of in_tree node
	double ret=0;
	inc(i,2) 
		if (t.p[i]>s[u].y[i]) ret+=t.p[i]-s[u].y[i];
		else if (t.p[i]<s[u].x[i]) ret+=s[u].x[i]-t.p[i];
	return ret;
}
double disL1Max(int u, Node &t){
	double ret=0;
	inc(i,2) ret+=max(abs(t.p[i]-s[u].x[i]),abs(t.p[i]-s[u].y[i]));
	return ret;
}
double sqr(double a){
	return a*a;
}
double disL2Min(int u, Node &t){
	double ret=0;
	inc(i,2) 
		if (t.p[i]>s[u].y[i]) ret+=sqr(t.p[i]-s[u].y[i]);
		else if (t.p[i]<s[u].x[i]) ret+=sqr(t.p[i]-s[u].x[i]);
	return ret;
}
double disL2Max(int u, Node &t){ //max coord dis
	double ret=0;
	inc(i,2) ret+=max(sqr(t.p[i]-s[u].x[i]),sqr(t.p[i]-s[u].y[i]));
	return ret;
}
void build(int &u, int l, int r, int cur){ //O(nlogn)
	u=l+r>>1; K=cur;
	nth_element(a+l,a+u,a+r+1);
	add(u,a[u]);
	if (l<u) build(lc,l,u-1,cur^1);
	if (r>u) build(rc,u+1,r,cur^1);
	upd(u);
}
//Maybe we need to rebuild the tree after unbalanced insert
void ins(int u, int cur){  
	if (q.p[cur]<s[u].p[cur])
		if (lc) ins(lc,cur^1);
		else lc=++n,add(n,q);
	else
		if (rc) ins(rc,cur^1);
		else rc=++n,add(n,q);
	upd(u);
}
void ask(int u){
	double di=sqr(s[u].p[0]-q.p[0])+sqr(s[u].p[1]-q.p[1]);
	if (di<ans[1]){
		ans[1]=di;
		ansp[1]=s[u].id;
	}
	else if (di<ans[0]){
		ans[0]=di;
		ansp[0]=s[u].id;
	}
	double dl=lc?disL2Min(lc,q):inf, dr=rc?disL2Min(rc,q):inf;
	//int dl=lc?disL1Max(lc,q):0, dr=rc?disL1Max(rc,q):0;
	if (dl<dr){ //trim branch, swap > < when search max dis point
		if (dl<ans[0]) ask(lc);
		if (dr<ans[0]) ask(rc);
	}
	else{
		if (dr<ans[0]) ask(rc);
		if (dl<ans[0]) ask(lc);
	}
}
//minDisPoint (L2 dis) with ins operate
//each query asks one nearest point of a giving coord
int main(){
	freopen("CA.txt","r",stdin);
	while(~scanf("%*d%*d%*d%*f%*f%*f%lf%lf",&a[n].p[0],&a[n].p[1])) a[n].id=n, n++;
	clock_t t0=clock();
	printf("%d data point loaded in %lu ms\n",n, t0);
	printf("building k-d tree\n");
	build(rt,0,n-1,0);
	printf("k-d tree builded in %lu ms\n",clock()-t0);
	t0=clock();
	
	printf("trying %d questions\n",n);
	for (int i=0;i<n;i++){
		q=a[i];
		ans[0]=ans[1]=inf; ask(rt);
		//printf("%d %d %f\n",a[i].id,ansp[0],ans[0]);
	}
	printf("done %d questions in %lu ms\n",n,clock()-t0);
	return 0;
}

