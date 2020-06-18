/*
Author: fffasttime
Date: 2020/06/12
Description: Assiociation rules discovery on US house voting dataset, using Apriori algorithm 
*/

#include <iostream>
#include <sstream>
#include <fstream>
#include <map>
#include <vector>
#include <string>
#include <bitset>
#include <algorithm>
using namespace std;

typedef unsigned long long ll;
const int N=435, M=34;

string AttrName[]={
"democrat",
"handicapped-infants",
"water-project-cost-sharing",
"adoption-of-the-budget-resolution",
"physician-fee-freeze",
"el-salvador-aid",
"religious-groups-in-schools",
"anti-satellite-test-ban",
"aid-to-nicaraguan-contras",
"mx-missile",
"immigration",
"synfuels-corporation-cutback",
"education-spending",
"superfund-right-to-su",
"crime",
"duty-free-exports",
"export-administration-act-south-africa"
};

ll data[N];

map<ll, int> fset;

const int Freq=150;
const double Confidence=0.8;

vector<pair<ll,ll>> rules;

void prt(ll x){
	cout<<(bitset<M>(x));
}

double lift(ll f, ll u, ll v){
	return (double)N*fset[f]/fset[u]/fset[v];
}

// u union v = f, u intersect v = empty
void dfs(ll f, ll u){
	ll v=f^u;
	if ((double)fset[f]/fset[v]>Confidence){
		if (lift(f,u,v)>2.0)
			rules.emplace_back(f, v);
	}
	else 
		return;
	for (int i=0;i<M;i++)
		if ((~u>>i&1) && (f>>i&1))
			dfs(f,u|1ull<<i);
}

void prtattr(ll x){
	for (int i=0;i<M/2;i++)
		if (x>>(i*2)&1)
			cout<<AttrName[i]<<":y ";
		else if (x>>(i*2+1)&1)
			cout<<AttrName[i]<<":n ";
}

int count(ll x){
	int cnt=0;
	for (int i=0;i<N;i++)
		if ((x&data[i]) == x)
			cnt++;
	return cnt;
}

int main(){
	ifstream fin("house-votes-84.data");
	for (int i=0;i<N;i++){
		string line;
		getline(fin,line);
		istringstream is(line);
		for (int j=0;j<M/2;j++){
			string s;
			getline(is,s,',');
			if (s=="democrat" || s=="y")
				data[i]|=1ull<<(j*2);
			else if (s!="?")
				data[i]|=1ull<<(j*2+1);
		}
	}
	map<ll, int> cur;
	cur[0]=fset[0]=N;
	for (int i=0;i<M;i++){ //get frequecny set
		map<ll, int> next;
		for (auto p:cur){
			ll x=p.first;
			for (int i=0;i<M;i++) if (~x>>i&1){
				int c=count(x|1<<i);
				if (c>Freq)
					next[x|1<<i]=c;
			}
		}
		if (next.empty()) break;
		cur=next;
		fset.insert(begin(next),end(next));
	}
	for (auto f:fset){ 
		cout<<f.second<<' ';
		prt(f.first); cout<<'\n';
	}
	for (auto f:fset){ //generate rules
		dfs(f.first,0);
	}
	sort(begin(rules),end(rules));
	rules.erase(unique(begin(rules),end(rules)),end(rules));
	cout<<fset.size()<<" frequent items\n"<<"Found "<<rules.size()<<" rules"<<"\n";
	for (auto p:rules){
		prt(p.second);
		cout<<"->";
		prt(p.first);
		cout<<" "<<(double)fset[p.first]/fset[p.second]<<'\n';
		prtattr(p.second);
		cout<<"=>\n";
		prtattr(p.first);
		cout<<"lift:"<<lift(p.first, p.second, p.first^p.second)<<"\n\n";
	}
	return 0;
}

