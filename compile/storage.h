/**
 * This header file defines storage method of data
*/

#ifndef STORAGE_H
#define STORAGE_H

#include "common.h"
#include <cassert>
#include <cstdio>
#include <cstring>
#include <vector>
#include <map>
using std::vector;
using std::map;
using std::FILE;

class IOException:exception{};

struct FieldReturn_t{
    void *p;
    int &int32(){return *(int *)p;}
    char &int8(){return *(char *)p;}
    char *nchar(){return (char *)p;}
    float &float32(){return *(float *)p;}
    //set a string
    void nchar(char *s){
        strcpy((char *)p, s);
    }
};

struct FieldCellInfo{
    FieldType type;
    int extra; // some type have extra info
    string name;
    char constraint; // 1 NOT NULL | 2 UNIQUE
    bool allowNULL() const {return constraint==0;}
    int length() const{
        switch (type){
        case FieldType::int32:
            return 4;
        case FieldType::nchar:
            // '\0' need one more byte to save
            return extra + 1;
        case FieldType::int8:
            return 1;
        case FieldType::Float:
            return 4;
        default:
            fatal("fatal: error loading field type %d, the file is broken?\n", int(type));
        }
    }
    FieldCellInfo(FieldType ft, string name):type(ft),name(name){}
    FieldCellInfo(FieldType ft, int extra, string name):type(ft),extra(extra),name(name){}
    FieldCellInfo(FieldType ft, int extra, string name, char constraint):
        type(ft),extra(extra),name(name),constraint(constraint){}
};

struct FieldInfo{
    vector<int> offset; //for faster calcuate pos
    map<string, int> name2fid; //find by name
    vector<FieldCellInfo> fields;
    int length;

    FieldInfo(const vector<FieldCellInfo> &__v);
};

typedef char* PRecord_t;

class Table{
public:
    FieldInfo field;
    bool loaded;
    vector<PRecord_t> data;
    string name;
    
    void loadData(FILE *fi, string _name, int _item_count);
    void saveData(FILE *fo);
    void freeData();
    void showData();
    
    FieldReturn_t read(int line, string fieldname){
        #ifndef NODEBUG
            assert(field.name2fid.count(fieldname));
        #endif
        return {data[line]+field.offset[field.name2fid[fieldname]]};
    }
    FieldReturn_t readof(int line, int offset){
        #ifndef NODEBUG
            assert(offset>=0 && offset<1<<16);
        #endif
        return {data[line]+offset};
    }
    int getoffset(string fieldname){
        return field.offset[field.name2fid[fieldname]];
    }

    Table(const FieldInfo &_field):field(_field),loaded(false){}
    ~Table();
};

/* Read-only version*/
class TableSt{
    FieldInfo field;
    vector<PRecord_t> data;
    PRecord_t data0;
    TableSt(const FieldInfo &_field, int lines):field(_field){
        data0 = (PRecord_t)malloc(lines*field.length);
        data.reserve(lines);
        for (int i=0;i<lines;i++)
            data.push_back(data0+i*field.length);
    }
    ~TableSt(){
        free(data0);
    }
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
