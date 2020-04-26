/*
Author: fffasttime
Date: 2020/04/23
Description: A simple demo of particle swarm optimization(PSO) algorithm
*/

#include <bits/stdc++.h>
using namespace std;

typedef double db;

const db PI=acos(-1);
struct Point{
	db x,y;
	Point operator+(const Point &v)const{return {x+v.x,y+v.y};}
	Point operator-(const Point &v)const{return {x-v.x,y-v.y};}
	Point operator*(db a)const{return {a*x,a*y};}
};

mt19937 rand_gen(0);

// accurate result: x* =11.6255, y* =5.7250, f(x*, y*)=38.8503
const db XMIN=-3.0, XMAX=12.1; 
const db YMIN=4.1, YMAX=5.8; //bound
db f(Point p){
	db x=p.x, y=p.y;
	return (21.5 + x*sin(4*PI*x) + y*sin(20*PI*y))-(x<XMIN || x>XMAX || y<YMIN || y>YMAX)*10000;
}

//update formula: p(k+1)=w*p(k) + c1 * pbest + c2 * gbest
const db WBEGIN=0.9, WEND=0.4, C1=1, C2=1;
const int ITER=500;

const int POP_SIZE=50;
struct P{
	Point p, pbest, v;
	db fitness;
	P(){
		uniform_real_distribution<db> disx(XMIN,XMAX);
		uniform_real_distribution<db> disy(YMIN,YMAX);
		uniform_real_distribution<db> disv(-1,1);
		p.x=disx(rand_gen); p.y=disy(rand_gen);
		v.x=disv(rand_gen); v.y=disy(rand_gen);
		pbest=p;
		fitness=f(p);
	}
	void upd(db w, Point gbest){
		uniform_real_distribution<db> dis01(0,1);
		v=v*w + (pbest-p)*C1*dis01(rand_gen) + (gbest-p)*C2*dis01(rand_gen);
		p=p+v;
		if (f(p)>f(pbest)) pbest=p;
		fitness=f(p);
	}
}g[POP_SIZE];

void update(int iter){
	db w = WBEGIN + (WEND - WBEGIN) * (1.0*iter/ITER); //linear
	Point gbest = max_element(g,g+POP_SIZE,[](const P& a, const P &b){return a.fitness<b.fitness;})->pbest;
	for (int i=0;i<POP_SIZE;i++){
		g[i].upd(w, gbest);
	}
}

void debug(int iter){
	if (iter%100) return;
	printf("iter %d:\n",iter);
	for (int i=0;i<POP_SIZE;i++){
		printf("(%.5f, %.5f), %.5f\n", g[i].p.x, g[i].p.y, g[i].fitness);
	}
}

int main(){
	debug(0);
	for (int i=1;i<=ITER;i++){
		update(i);
		debug(i);
	}
	return 0;
}

