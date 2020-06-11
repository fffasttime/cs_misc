/*
Author: fffasttime
Date: 2020/05/26 
Description: LSH (Locality sensitive hashing), an algorithm 
for approximate nearest neighbors search problem
*/
#include <set>
#include <map>
#include <algorithm>
#include <ctime>
#include <random>
#include "dataload.hpp"
using namespace std;

DataVec randvec(){
	DataVec r;
	for (int i=0;i<D;i++)
		r[i]=rand()-RAND_MAX/2;
	db l=L2(r);
	for (int i=0;i<D;i++)
		r[i]/=l;
	return r;
}

Data data;
const int Round=5,Bit=12;
	
struct Hash{
	multimap<int,int> h;
	vector<DataVec> dir; //direction vector
	int gethash(const DataVec &d){
		int id=0;
		for (int j=0;j<Bit;j++){
			db x=dot(dir[j], d);
			//TODO: add constant item b
			id<<=1;
			id|=x>0;
		}
		return id;
	}
	void construct(){
		dir.resize(Bit);
		for (int i=0;i<Bit;i++) dir[i]=randvec();
		for (int i=0;i<N;i++)
			h.insert({gethash(data[i]),i});
	}
	vector<int> search(const DataVec &v){
		auto range=h.equal_range(gethash(v));
		vector<int> ret;
		for (auto it=range.first;it!=range.second;it++)
			ret.push_back(it->second);
		return ret;
	}
};

int main()
{
	Hash h[Round]; // hash table
	clock_t t0=clock();
	for (int i=0;i<Round;i++)
		h[i].construct();
	cout<<"constructed hash table\n";
	clock_t t1=clock();
	
	int acc=0,mcc=0;
	//test
	int nearsum=0;
	for (int i=0;i<cN;i++){
		set<int> nears;
		for (int j=0;j<Round;j++){
			auto ret=h[j].search(data[i]);
			nears.insert(begin(ret),end(ret));
		}
		vector<pair<db,int>> near;
		for (auto v: nears) near.emplace_back(L2(data[i],data[v]),v);
		check(nears.size()>10, "near set too small!");
		nearsum+=near.size();
		partial_sort(begin(near),begin(near)+11,end(near));
		
		bool matched=1;
		for (int j=1;j<min(11,(int)near.size());j++){
			if (near[j].second!=data.real[i][j-1])
				matched=0;
			else mcc++;
		}
		acc+=matched;
	}
	cout<<"average hash near set:"<<(db)nearsum/cN<<'\n';
	FILE *res=fopen("result.md","a");
	fprintf(res,"|r=%2d|len=%3d|acc=%5.2f%%|recall=%5.2f%%|%4ld ms|%5ld ms|\n",
		Round, Bit, (db)acc/cN*100, (db)mcc/cN*10, t1-t0, clock()-t1);
	fclose(res);
	
	return 0;
}

