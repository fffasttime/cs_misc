/*
Author: fffasttime
Date: 2020/03/15
*/

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <cstring>
#include <fstream>
#include <vector>
#include <iostream>
using namespace std;

void check(bool pred, const char *reason){
	if (!pred){
		printf("%s\n", reason);
		exit(1);
	}
}

string filename;

double range[129];

void statRange(char *content, int len){
	int count[128];
	memset(count,0,sizeof count);
	for (int i=0;i<len;i++) 
		count[content[i]]++;
	//for (int i=0;i<128;i++) if (count[i]==1) count[i]++,len++;
	FILE *dict=fopen((filename+".dict").c_str(),"w");
	for (int i=0;i<128;i++) fprintf(dict,"%d\n",count[i]);
	//for (int i=0;i<128;i++) if (count[i]) cout<<i<<' '<<count[i]<<'\n';
	fclose(dict);
	for (int i=0;i<128;i++) count[i+1]+=count[i];
	for (int i=0;i<128;i++) range[i+1]=(double)count[i]/count[127];
	for (int i=0;i<128;i++) cout<<range[i]<<'\n';
}

typedef unsigned long long uint;
typedef unsigned char byte;
const uint P=32; //upper range 

void encode(){
	FILE *input=fopen(filename.c_str(),"r");
	check(input, "can't find source file");
	auto content=new char[1000000];
	content[fread(content,1,1000000,input)]=0;
	fclose(input);
	int len=strlen(content)+10; //padding
	statRange(content,len);
	
	//initial
    uint L = 0, R = 1ull << P;
    byte buf = -1, c = 0, z = 0; bool first=1;
    ofstream out(filename+".bin", ios::binary);
	
	for (int i=0;i<len;i++){
		int ch=content[i];
		//update
	    L = L + R * (double)range[ch];
	    R = R * (double)(range[ch+1] - range[ch]);
	    uint temp = L >> (P - 8);
	    if (temp > 0xff){ //overflow?
	        buf++ ;
	        L -= (1ull << P);
	        z = 1;
	    }
	    
	    //renormalize
	    while(R < (1 << (P - 8))){
			temp = L >> (P - 8);
			if (temp < 0xff){
				if (!first) out.put(buf);   //output previous 
				for (;c > 0;c--)   // output 0xff or 0x00
					if(z == 1)
						out.put(0);
					else
						out.put(0xff);
				buf = temp; // current buffer
				first=0;
			}
			else if (temp == 0xff) // can't judge now
				c++;
			R <<= 8;
			L = (L << 8) & ((1ull << P) - 1); // left shift , and cut higher positon
			z = 0;
    	}
	}
	//out.put(L >> (P - 8)); //last
	delete[] content;
	cout<<"Encode done in "<<filename+".bin"<<endl;
}
void decode(){
	FILE *dict=fopen((filename+".dict").c_str(),"r");
	check(dict, "can't find dict file");
    ifstream in(filename+".bin", ios::binary);
	check(in.is_open(), "can't find binary file");
	
	int count[128];
	for (int i=0;i<128;i++) fscanf(dict, "%d", count+i);
	fclose(dict);
	for (int i=0;i<128;i++) count[i+1]+=count[i];
	for (int i=0;i<128;i++) range[i+1]=(double)count[i]/count[127];
	
	uint L = 0, R = 1ull << P, V=in.get();
	V=V<<8 | in.get();
	V=V<<8 | in.get();
	V=V<<8 | in.get();
    
    ofstream out(string("out_")+filename);
    
	while (!in.eof()){
		// get symbol
		if(V < L) // carry bit 
			V += 1ull << P;
    	double T = (double)(V - L) / R;
		int ch;
		for (int i=0;i<128;i++) 
			if(range[i+1]-range[i]>=1e-8 && T>=range[i] && T<range[i+1]){
				ch=i; break;
			}
		if (ch==0) break;
		// update
		out<<(char)ch;
	    L = L + R * (double)range[ch];
	    //cout<<R;
	    //cout<<' '<<R * (double)(range[ch+1] - range[ch])<<'\n';
	    R = R * (double)(range[ch+1] - range[ch]);
	    while(R < (1 << (P - 8))){
        	R <<= 8;
        	L = (L << 8) & ((1ull << P) - 1);
        	V = (V << 8) & ((1ull << P) - 1);
    		V |= in.get();
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
