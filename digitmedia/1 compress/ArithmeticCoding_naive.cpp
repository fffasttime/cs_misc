/*
Author: fffasttime
Date: 
*/

// only store in one double variable
// so input lenth should less then 4
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <cstring>
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
	vector<pair<int,double>> prob;
	for (int i=0;i<128;i++){
		range[i+1]=range[i]+(double)count[i]/len;
		if (count[i]) cout<<i<<' '<<count[i]<<'\n';
	}
	FILE *dict=fopen((filename+".dict").c_str(),"w");
	for (int i=0;i<129;i++)
		fprintf(dict,"%.15f\n",range[i]);
	fclose(dict);
}

void zip(){
	FILE *input=fopen(filename.c_str(),"r");
	check(input, "can't find source file");
	auto content=new char[1000000];
	fread(content,1,1000000,input);
	fclose(input);
	int len=strlen(content);
	statRange(content,len);
	double low=0.0,up=1.0;
	for (int i=0;i<len;i++){
		int ch=content[i];
		double section=up-low;
		up=low+range[ch+1]*section;
		low=low+range[ch]*section;
	}
	printf("%.15f\n",low);
	delete[] content;
	
}
void unzip(){
	FILE *dict=fopen((filename+".dict").c_str(),"r");
	check(dict, "can't find dict file");
	for (int i=0;i<129;i++)
		fscanf(dict, "%lf",range+i);
	fclose(dict);
	double x; cin>>x;
	while (x>0){
		for (int i=0;i<128;i++)
			if (x>=range[i] && x<range[i+1]){
				x-=range[i];
				x/=range[i+1]-range[i];
				cout<<(char)i;
				break;
			}
	}
}

int main(int argc, char **argv){
	check(argc==3, "bad args");
	filename=string(argv[2]);
	if (argv[1][0]=='u') unzip();
	else if (argv[1][0]=='z') zip();
	else check(0,"unzip or zip?");
	return 0;
}
