/**
 * This header file defines storage method of data
*/

#ifndef STORAGE_H
#define STORAGE_H

#include "common.h"
#include <cstdio>
#include <cstring>
#include <vector>
#include <map>
using std::vector;
using std::map;
using std::FILE;

struct FieldReturn_t{
    void *p;
    int &int32(){return *(int *)p;}
    char *nchar(){return (char *)p;}
    //set a string
    void nchar(char *s){
        strcpy((char *)p, s);
    }
};

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
            fatal("fatal: error loading field type\n");
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

    FieldInfo(const vector<FieldCellInfo> &__v):fields(__v){
        offset.resize(fields.size());
        for (size_t i=1;i<offset.size();i++)
            offset[i]=offset[i-1]+fields[i-1].length();
        length=offset.back()+fields.back().length();
        int fid=0;
        for (const auto &v:fields)
            name2fid[v.name]=fid++;
    }
};

typedef char* PRecord_t;

class Table{
public:
    int count;
    FieldInfo field;
    bool loaded;
    vector<PRecord_t> data;
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

    Table(const FieldInfo &_field):field(_field),loaded(false){}
    ~Table();
};

FieldInfo getDBMetaFieldInfo();
FieldInfo getTableMetaFieldInfo();

class DataBase{
public:
    string path;
    string name;
    bool loaded;
    Table dbmeta;
    vector<Table> tables;
    vector<Table> tbmeta;
    //update when tables change
    map<string, Table*> name_tab;
private:
    void freeData();
    void loadData_tables(FILE *fi, int table_count);
    bool checkDBFile(FILE *fi, int &table_count);
public:
    void updateMap();
    DataBase():loaded(false),dbmeta(getDBMetaFieldInfo()){}
    ~DataBase();
    void loadData(string _path, string _name);
    void saveData();
    void createNew(string _path, string _name);
};

#endif
