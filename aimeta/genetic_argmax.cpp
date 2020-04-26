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
const int N = 50; //population size 
const int E=2; //number of new_gr alive
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

int random_choice(const vector<db> &c){
	uniform_real_distribution<db> dis(0,c.back());
	return lower_bound(begin(c),end(c),dis(rand_gen))-begin(c);
}

struct Gene{
	ll code;
	Gene(){}
	Gene(ll c):code(c){}
	pair<db,db> decode() const{
		db p1=code>>GENE_SECOND_PART;
		db p2=code&((1ll<<GENE_SECOND_PART)-1);
		return {XMIN + (XMAX-XMIN) * p1/(1ll<<GENE_FIRST_PART), 
				YMIN + (YMAX-YMIN) * p2/(1ll<<GENE_SECOND_PART)};
	}
	db fitness() const{
		auto p=decode();
		//if (p.first<XMIN || p.first>XMAX){printf("%f\n",p.first);cout<<bitset<35>(code);}
		return target(p.first, p.second);
	}
};
vector<Gene> gene;

void init(){
	gene.clear();
	uniform_int_distribution<ll> dis(0,(1ll<<GENE_LEN)-1); //close section
	for (int i=0;i<N;i++) gene.emplace_back(dis(rand_gen));
}

void iter(){
	vector<db> v,c;
	for (int i=0;i<N;i++) v.push_back(gene[i].fitness());
	partial_sum(begin(v),end(v),back_inserter(c));
	vector<int> cross_pool;
	for (int i=0;i<N;i++)
		cross_pool.push_back(random_choice(c));
	uniform_real_distribution<db> dis01(0,1);
	vector<Gene> new_gr;
	//cross operator
	for (int i=0;i<N;i+=2)
		if (dis01(rand_gen)<P_C){
			ll c1=gene[cross_pool[i]].code;
			ll c2=gene[cross_pool[i+1]].code;
			uniform_int_distribution<> dis(0,GENE_LEN-1);
			int l=dis(rand_gen), r=dis(rand_gen);
			bool ff=0;
			if (l>r) swap(l,r),ff=1;
			ll bitmask=((1ull<<r)-1)>>l<<l; //section [l,r)
			//if (ff) bitmask^=(1ll<<GENE_LEN)-1;
			new_gr.emplace_back(c1&bitmask | c2&~bitmask);
			new_gr.emplace_back(c2&bitmask | c1&~bitmask);
		}
	//mutation operator
	int nsize=new_gr.size();
	for (int i=0;i<nsize;i++){
		ll c0=new_gr[i].code;
		ll c1=c0;
		for (int j=0;j<GENE_LEN;j++)
			if (dis01(rand_gen)<P_M)
				c1^=1ll<<j;
		if (c1!=c0)
			new_gr.emplace_back(c1);
	}
	vector<Gene> next_gr;
	if (new_gr.size()<E) new_gr.push_back(gene[0]);
	if (new_gr.size()<E) new_gr.push_back(gene[1]);
	//for (int i=0;i<N;i++) new_gr.push_back(gene[i]);
	partial_sort(begin(new_gr),begin(new_gr)+E,end(new_gr),
		[](const Gene&x, const Gene&y){return x.fitness()>y.fitness();});
	for (int i=0;i<N;i++)
		if (i<E) next_gr.push_back(new_gr[i]);
		else next_gr.push_back(gene[random_choice(c)]);
	
	swap(gene, next_gr);
}

void debug(int iter){
	if (iter%500==0){
		printf("iteration %d:\n",iter);
		for (int i=0;i<N;i++){
			auto p=gene[i].decode();
			printf("(%2.5f, %2.5f), fitness=%f\n", p.first, p.second, gene[i].fitness());
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
			iter();
			//debug(i);
		}
		if (gene[0].fitness()>38.85) auc++;
	}
	printf("|%.2fms|%.2f%% (%d/%d)|\n",clock()*1.0/CT, auc*100.0/CT, auc, CT);
	return 0;
}

