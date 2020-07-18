/**
 * This header file defines storage method of data
*/

#ifndef STORAGE_H
#define STORAGE_H

#include "util.h"
#include <stdio.h>
#include <vector>
#include <map>
using std::vector;
using std::map;

struct FieldReturn_t{
    void *p;
    int int32(){return *(int *)p;}
    char *nchar(){return *(char *)p;}
}

struct FieldCellInfo{
    FieldType type;
    int extra; //some type have extra info
    string name;
    int length(){
        switch (type){
        case FieldType::int32:
            return 4;
        case FieldType::nchar:
            return extra;
        default:
            fatal("Fatal: error loading field type\n");
        }
    }
    FieldCellInfo(FieldType ft, string name):type(ft),name(name){}
    FieldCellInfo(FieldType ft, int extra, string name):type(ft),extra(extra),name(name){}
};

struct FieldInfo{
    vector<int> offset; //for faster calcuate pos
    map<string, int> name2fid; //find by name
    vector<FieldCellInfo> fields;
    int length;

    FieldInfo(vector<FieldCellInfo> __v):fields(__v){
        offset.resize(fields.size());
        for (int i=1;i<offset.size();i++)
            offset[i]=offset[i-1]+fields[i-1].length();
        int fid=0;
        for (const auto &v:fields)
            name2fid[v.name]=fid++;
    }
}

typedef char* PRecord_t;

class Table{
public:
    int tid;
    int item_count;
    bool loaded;
    vector<PRecord> data;
    FieldInfo field;
    string name;
    
    void loadData(FILE *fi, string _name, int _item_count);
    void saveData(FILE *fo);
    void freeData();

    FieldReturn_t read(int line, string fieldname){
        return {data[line]+field.offset[field.name2fid[fieldname]]};
    }
    FieldReturn_t readof(int line, int offset){
        return {data[line]+offset};
    }
    int getoffset(string fieldname){
        return field.offset[field.name2fid[fieldname]];
    }

    Table(_field):field(_field),loaded(false){}
    ~Table();
};

class DataBase{
private:
    string path;
    string name;
    bool loaded;
    vector<Table> tables;
    map<string> name2id()
private:
    freeData();
public:
    DataBase():loaded(false){}
    ~DataBase();
    loadData(string _path, string _name);
    saveData();
};

#endif