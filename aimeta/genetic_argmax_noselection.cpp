/*
Author: fffasttime
Date: 2020/04/20
Description: A simple demo of genetic algorithm 
*/
#include <bits/stdc++.h>
using namespace std;

mt19937 rand_gen(0);
typedef unsigned long long ll;
typedef double db;

const db PI = acos(-1);
const int N = 50, E=30; //population size 
const int GENE_LEN = 33, GENE_FIRST_PART = 18, GENE_SECOND_PART = GENE_LEN - GENE_FIRST_PART;
const db P_C = 0.6; //cross prob
const db P_M = 0.1; //mutation prob

// accurate result: x* =11.6255, y* =5.7250, f(x*, y*)=38.8503
// it's a positive function, so can be used as fitness function directly 
const db XMIN=-3.0, XMAX=12.1; 
const db YMIN=4.1, YMAX=5.8; //bound
db target(db x, db y){
	assert(x>=XMIN && x<=XMAX);
	assert(y>=YMIN && y<=YMAX);
	return 21.5 + x*sin(4*PI*x) + y*sin(20*PI*y);
}

struct Gene{
	ll code;
	db fitness;
	Gene(){}
	Gene(ll c):code(c){calcfitness();}
	pair<db,db> decode() const{
		db p1=code>>GENE_SECOND_PART;
		db p2=code&((1ll<<GENE_SECOND_PART)-1);
		return {XMIN + (XMAX-XMIN) * p1/(1ll<<GENE_FIRST_PART), 
				YMIN + (YMAX-YMIN) * p2/(1ll<<GENE_SECOND_PART)};
	}
	void calcfitness(){
		auto p=decode();
		//if (p.first<XMIN || p.first>XMAX){printf("%f\n",p.first);cout<<bitset<35>(code);}
		fitness=target(p.first, p.second);
	}
};
vector<Gene> gene;

void init(){
	gene.clear();
	uniform_int_distribution<ll> dis(0,(1ll<<GENE_LEN)-1); //close section
	for (int i=0;i<N;i++) gene.emplace_back(dis(rand_gen));
}

void iter(int it){
	//cross operator
	uniform_real_distribution<db> dis01(0,1);
	uniform_int_distribution<> disN(0,N-1);
	for (int i=0;i<N;i+=2)
		if (dis01(rand_gen)<P_C){
			ll c1=gene[disN(rand_gen)].code;
			ll c2=gene[disN(rand_gen)].code;
			uniform_int_distribution<ll> disb(0,(1ll<<GENE_LEN)-1); //close section
			ll bitmask=disb(rand_gen);
			gene.emplace_back(c1&bitmask | c2&~bitmask);
			gene.emplace_back(c2&bitmask | c1&~bitmask);
		}
	//mutation operator
	int nsize=gene.size();
	for (int i=N;i<nsize;i++){
		ll c0=gene[i].code;
		ll c1=c0;
		for (int j=0;j<GENE_SECOND_PART;j++)
			if (dis01(rand_gen)<P_M)
				c1^=1ll<<j;
		for (int j=GENE_SECOND_PART;j<GENE_LEN;j++)
			if (dis01(rand_gen)<P_M)
				c1^=1ll<<j;
		if (c1!=c0)
			gene.emplace_back(c1);
	}
	nth_element(begin(gene),begin(gene)+E,end(gene),
		[](const Gene&x, const Gene&y){return x.fitness>y.fitness;});
	shuffle(begin(gene)+E,end(gene),rand_gen);
	gene.resize(N);
}

void debug(int iter){
	if (iter%500==0){
		printf("iteration %d:\n",iter);
		for (int i=0;i<N;i++){
			auto p=gene[i].decode();
			printf("(%2.5f, %2.5f), fitness=%f\n", p.first, p.second, gene[i].fitness);
		}
	}
}

int main(){
	int CT=500, auc=0;
	for (int T=0;T<CT;T++){
		init();
		const int ITER=500;
			//debug(0);
		for (int i=1;i<=ITER;i++){
			iter(i);
			//debug(i);
		}
		if (gene[0].fitness>38.85) auc++;
	}
	printf("|%.2fms|%.2f%% (%d/%d)|\n",clock()*1.0/CT, auc*100.0/CT, auc, CT);
	return 0;
}

