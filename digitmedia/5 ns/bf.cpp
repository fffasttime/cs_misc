#include <iostream>
#include <cstdio>
#include <cmath>
#include <cstring>
#include <ctime>
#include <algorithm>
using namespace std;

const double inf=1e10;
double ans; int ansp;

int K,n;
struct Node{
	double p[2],x[2],y[2]; int id;
	bool operator<(const Node &v)const{
		return p[K]<v.p[K];
	}
}a[100000];

double sqr(double x){
	return x*x;
}

double l2dis(int x, int y){
	return sqr(a[x].p[0]-a[y].p[0])+sqr(a[x].p[1]-a[y].p[1]);
}

void ask(int u){
	for (int i=0;i<n;i++)
		if (i-u && ans>l2dis(i,u)){
			ans=l2dis(i,u);
			ansp=i;
		}
}

int main(){
	FILE *in=fopen("CA.txt","r");
	while(~fscanf(in,"%*d%*d%*d%*f%*f%*f%lf%lf",&a[n].p[0],&a[n].p[1])) a[n].id=n, n++;
	clock_t t0=clock();
	printf("%d data point loaded in %lu ms\n",n, t0);
	//int x,y; while (1){cin>>x>>y;cout<<l2dis(x,y)<<'\n';}
	
	t0=clock();
	
	printf("trying %d questions\n",n);
	for (int i=0;i<n;i++){
		ans=inf; ask(i);
		//printf("%d %d %f\n",a[i].id, ansp, sqrt(ans));
	}
	printf("done %d questions in %lu ms\n",n,clock()-t0);
	return 0;
}
