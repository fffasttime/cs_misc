/*
Author: fffasttime
Date: 2020/03/19
*/

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <cstring>
#include <fstream>
#include <vector>
#include <map>
#include <iostream>
using namespace std;

void check(bool pred, const char *reason){
	if (!pred){
		printf("%s\n", reason);
		exit(1);
	}
}

string filename;


typedef unsigned long long uint;
typedef unsigned char byte;
const uint P=32; //upper range 

void encode(){
	FILE *input=fopen(filename.c_str(),"r");
	check(input, "can't find source file");
	auto content=new char[1000000];
	content[fread(content,1,1000000,input)]=0;
	fclose(input);
	int len=strlen(content); //padding
	
	ofstream out(filename+".bin", ios::binary);
    //vector<string> dix;
    int clen=0;
    map<string,int> dict; int dicc=0;
    for (;dicc<128;dicc++) dict[string({(char)dicc})]=dicc;
    string P;
	for (int i=0;i<=len;i++){
    	char C=content[i];
    	if (!dict.count(P+C)){
    		//dix.push_back(P+C);
			dict[P+C]=dicc++;
    		int x=dict[P];
    		if (x<128) out.put(x),clen++;
    		else if(x<16384){
    			out.put(x>>8|0x80);
    			out.put(x&0xFF);
    			clen+=2;
			}
			else{
    			out.put(x>>16|0xC0);
    			out.put(x>>8&0xFF);
    			out.put(x&0xFF);
    			clen+=3;
			}
			P=C;
		}
		else P+=C;
	}
    
    //{ofstream out(filename+".dict");for (auto &s: dix) out<<s<<'\n';}
    
	cout<<"LZW encoding lenth: "<<clen<<" byte"<<endl;
	delete[] content;
	cout<<"Encode done in "<<filename+".bin"<<endl;
}
void decode(){
    ifstream in(filename+".bin", ios::binary);
	check(in.is_open(), "can't find binary file");
	
    
    ofstream out(string("out_")+filename);
	vector<string> dict;
	for (int i=0;i<128;i++) dict.push_back(string({(char)i}));
    int cW=in.get();
    out<<dict[cW];
	while (1){
		int pW=cW;
		cW=in.get();
		if (in.eof()) break;
		if (cW&0x80) cW=(cW^0x80)<<8 | in.get();
		if (cW&0xC000) cW=(cW^0x4000)<<8 | in.get();
		if (cW<dict.size()){
			out<<dict[cW];
			dict.push_back(dict[pW]+dict[cW][0]);
		}
		else{
			string p=dict[pW];
			p=p+p[0];
			cW=dict.size();
			dict.push_back(p);
			out<<p;
		}
	}
	
	cout<<"Decode done in out_"<<filename<<endl;
}

int main(int argc, char **argv){
	check(argc==3, "bad args");
	filename=string(argv[2]);
	if (argv[1][0]=='d') decode();
	else if (argv[1][0]=='e') encode();
	else check(0,"encode or decode?");
	return 0;
}
