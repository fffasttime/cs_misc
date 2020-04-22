/*
Author: fffasttime
Date: 2020/04/23
Description: A simple demo of genetic algorithm solving 0-1 backpack problem
*/
#include <bits/stdc++.h>
using namespace std;

mt19937 rand_gen(0);
typedef unsigned long long ll;
typedef double db;

const db PI = acos(-1);
const int POP_SIZE = 100; //population size 
const int N = 50; // number of backpack inventories
const int GENE_LEN = N;
const db P_C = 0.6; //cross prob
const db P_M = 0.05; //mutation prob

//this test data: v* = 3103, sumw = 1000
int W0 = 1000;
db V[N]={220, 208,198,192,180,180,165,162,160,158,
155,130,125,122,120,118,115,110,105,101,
100,100,98,96,95,90,88,82,80,77,75,73,72,
70,69,66,65,63,60,58,56,50,30, 20,15,10,8,
5,3,1};
int W[N]={80,82,85,70,72,70,66,50,55, 25,50,55, 40, 48,
50,32, 22,60,30,32, 40,38,35,32, 25, 28,30, 22,
50,30, 45,30,60,50, 20,65, 20, 25,30,10, 20, 25,
15,10,10,10, 4, 4, 2,1};

int random_choice(const vector<db> &c){
	uniform_real_distribution<db> dis(0,c.back());
	return lower_bound(begin(c),end(c),dis(rand_gen))-begin(c);
}

struct Gene{
	ll code;
	db fitness;
	//Gene(){}
	Gene(ll c):code(c){ refine(); calcfit(); } //auto keeping available
	// make this genetic available by greedy dropping worthless inventory
	void refine(){
		int sumw=0;
		priority_queue<pair<db, int>> cur;
		for (int i=0;i<N;i++) if (code>>i&1){
			sumw+=W[i];
			cur.push({-V[i]/W[i], i});
		}
		while (sumw>W0){
			int out=cur.top().second;
			sumw-=W[out];
			code^=1ull<<out;
			cur.pop();
		}
	}
	db calcfit(){
		fitness=0;
		for (int i=0;i<N;i++) 
			if (code>>i&1)
				fitness+=V[i];
	}
};
vector<Gene> gene;

void init(){
	uniform_int_distribution<ll> dis(0,(1ll<<GENE_LEN)-1); //close section
	for (int i=0;i<POP_SIZE;i++) gene.emplace_back(dis(rand_gen));
	
	for (int i=0;i<N;i++){ //shuffle
		int id=rand()%N;
		//swap(V[i],V[id]); swap(W[i],W[id]);
	}
}

void iter(){
	vector<db> v,c;
	for (int i=0;i<POP_SIZE;i++) v.push_back(gene[i].fitness);
	partial_sum(begin(v),end(v),back_inserter(c));
	vector<int> cross_pool;
	for (int i=0;i<POP_SIZE;i++)
		cross_pool.push_back(random_choice(c));
	uniform_real_distribution<db> dis01(0,1);
	vector<Gene> new_gr;
	//cross operator
	for (int i=0;i<POP_SIZE;i+=2)
		if (dis01(rand_gen)<P_C){
			ll c1=gene[cross_pool[i]].code;
			ll c2=gene[cross_pool[i+1]].code;
			uniform_int_distribution<> dis(0,GENE_LEN-1);
			int l=dis(rand_gen), r=dis(rand_gen);
			if (l>r) swap(l,r);
			//ll bitmask=((1ull<<r)-1)>>l<<l; //section [l,r)
			uniform_int_distribution<ll> dis1(0,(1ll<<GENE_LEN)-1);
			ll bitmask=dis1(rand_gen);
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
	const int E=2; //number of new_gr alive
	if (new_gr.size()<E) new_gr.push_back(gene[0]);
	if (new_gr.size()<E) new_gr.push_back(gene[1]);
	partial_sort(begin(new_gr),begin(new_gr)+E,end(new_gr),
		[](const Gene&x, const Gene&y){return x.fitness>y.fitness;});
	for (int i=0;i<POP_SIZE;i++)
		if (i<E) next_gr.push_back(new_gr[i]);
		else next_gr.push_back(gene[random_choice(c)]);
	
	swap(gene, next_gr);
}

void debug(int iter){
	if (iter%100==0){
		printf("iteration %d:\n",iter);
		for (int i=0;i<POP_SIZE;i++){
			ll c=gene[i].code;
			ll sumw=0; for (int i=0;i<N;i++) if (c>>i&1) sumw+=W[i];
			printf("fitness=%f, sumw=%lld\n", gene[i].fitness, sumw);
		}
	}
}

int main(){
	init();
	const int ITER=500;
		debug(0);
	for (int i=1;i<=ITER;i++){
		iter();
		debug(i);
	}
	return 0;
}

