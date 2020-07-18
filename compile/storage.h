/**
 * This header file defines storage method of data
*/

#ifndef STORAGE_H
#define STORAGE_H

#include "util.h"
#include <stdio.h>

struct FieldInfo{
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
                fatal("Fatal: error loading field type");
            }
        }
    };
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

typedef void* PRecord;

class Table{
private:
    int tid;
    int item_count;
    bool loaded;
    vector<PRecord> data;
    FieldInfo field;

    void loaddata(FILE *fi);
    void savedata(FILE *fo);
    void freedata();
public:
    Table(/* args */);
    ~Table();
};

class DataBase{
private:
    string path;
    string dbname;
public:
    DataBase(string path, string dbname):path(path), dbname(dbname);
    ~DataBase();
    savedata();
}

#endif