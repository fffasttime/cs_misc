/*
Author: fffasttime
Date: 2020/04/20
Description: A simple demo of SA_algo, 
	easier than genetic algorithm on such a simple problem
*/

#include <bits/stdc++.h>
using namespace std;
const double PI=acos(-1);
struct Point{
	double x,y;
};
double E(Point p){
	double x=p.x, y=p.y;
	return -(21.5 + x*sin(4*PI*x) + y*sin(20*PI*y))+(x<-3 || x>12.1 || y<4.1 || y>5.8)*10000;
}
/*
argmin_SimulateAnneal(): get argmin E(x), while E is not a convex function.
To speed up, the energy function should be nearly smooth and partial convex.
Point: n-d vector, the dimension cannot be too large.
t0:start temprature
tend: end temprature
delta: temprature reduce factor
*/
double rand01() {return (rand()+0.5)/RAND_MAX;}
Point ans; double anse;
void argmin_SimulateAnneal(double t0=50, double tend=1e-5, double delta=0.999){
	Point p((Point){5,5});
	anse=1e12;
	double t=t0, ne=E(p); //current state energy
	while (t>=tend){
		Point p1((Point){p.x+(rand01()*t-t/2), p.y+(rand01()*t-t/2)});
		double te=E(p1), K=1;
		//cout<<te-ne<<' '<<ne<<' '<<t<<' '<<exp((ne-te)/t)<<' '<<ans.x<<'\n';
		if (te<ne || 0 && exp((ne-te)/t*K)>rand01()) //disabled jumpout
			p=p1, ne=te; //update
		if (ne<anse) ans=p, anse=ne; 
		//cout<<ans.x<<' '<<ans.y<<' '<<anse<<'\n';
		t*=delta;
	}
}
int main(){
	int CT=500, auc=0;
	srand(1);
	for (int T=0;T<CT;T++){
		argmin_SimulateAnneal();
		if (-anse>38.85) auc++;
		//printf("(%f %f): %f",ans.x,ans.y,-anse);
	}
	printf("|%.2fms|%.2f%% (%d/%d)|\n",clock()*1.0/CT, auc*100.0/CT, auc, CT);
}
