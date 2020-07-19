#include "storage.h"
#include <string.h>
#include <stdlib.h>

void Table::loadData(FILE *fi, string _name, int _item_count){
    assert(!loaded);
    data.resize(count);

    for (int i=0;i<count;i++){
        data[i] = (PRecord_t)malloc(field.length);
        fread(data[i], field.length, 1, fi);
    }

    loaded=1;
}
void Table::freeData(){
    assert(loaded);
    for (const auto &p : data)
        if (p != nullptr)
            free(p);
}

void Table::saveData(FILE *fo){
    assert(loaded);
    for (PRecord_t p : data)
        if (p!=nullptr){
            fwrite(p, field.length, 1, fo);
        }
}

Table::~Table(){
    if (loaded) 
        freeData();
}

/**
 * read first 64 bytes
 */
bool DataBase::checkDBFile(FILE *fi, int &table_count){
    char predata[64];
    if (fread(predata, 64, 1, fi) < 64){
        printf_error("error: bad db format of %s\n",path.c_str());
        return false;
    }
    // check head
    if (strncmp(predata,"simplesql",9)){
        printf_error("error: bad db format of %s\n",path.c_str());
        return false;
    }
    // check version
    const char cur_version[]="0.1";
    char file_version[10];
    strncpy(file_version, predata+10, 10);
    if (strncmp(file_version,"0.1",9)){
        printf_error("error: version check failed of %s, current engine is %s, but files %s\n", 
                        path.c_str(), cur_version, file_version);
        return false;
    }
    table_count=(int)predata[63];
    return true;
}

/**
 * Metadata table strucure, version 0.1
 */
FieldInfo getDBMetaFieldInfo(){
    vector<FieldCellInfo> fields;
    fields.emplace_back(FieldType::nchar,32,"name");
    fields.emplace_back(FieldType::int32,"count");
    fields.emplace_back(FieldType::int32,"count_field");
    return FieldInfo(fields);
}
FieldInfo getTableMetaFieldInfo(){
    vector<FieldCellInfo> fields;
    fields.emplace_back(FieldType::int32,"type");
    fields.emplace_back(FieldType::int32,"extra");
    fields.emplace_back(FieldType::nchar,32,"name");
    return FieldInfo(fields);
}

void DataBase::freeData(){
    loaded=false;
    tables.clear();
    name2tid.clear();
}

void DataBase::loadData_tables(FILE *fi, int table_count){
    tables.emplace_back(getDBMetaFieldInfo);
    Table &dbmeta=tables[0];
    dbmeta.loadData(fi, "__dbmeta", table_count);

    // load table metadata
    for (int i=0;i<table_count;i++){
        auto fieldinfo=getTableMetaFieldInfo();
        tables.emplace_back(fieldinfo);
        tables.back().loadData(fi, string("__tbmeta_")+dbmeta.read(i,"name").nchar(), 
                        dbmeta.read(i,"count_field").int32());
    }

    // load ordinary data
    for (int i=0;i<table_count;i++){
        vector<FieldCellInfo> fields;
        
        Table &tbmeta=tables[i+1];
        for (int j=0;j<tbmeta.count;j++)
        fields.emplace_back((FieldType)tbmeta.read(j,"type").int32(),
                            tbmeta.read(j,"extra").int32(), 
                            tbmeta.read(j,"name").nchar());
        tables.emplace_back(FieldInfo(fields));
        tables.back().loadData(fi, dbmeta.read(i,"name").nchar(), 
                        dbmeta.read(i,"count").int32());
    }
}

void DataBase::loadData(string _path, string _name){
    path=_path;
    name=_name;
    
    string fullpath=path+name+".db";
    
    FILE *fi=fopen(fullpath.c_str(), "rb");
    {
        if (!fi){
            printf_error("error: can't open file %s\n",path.c_str());
            return;
        }

        int table_count;
        if (!checkDBFile(fi, table_count)){
            fclose(fi);
            return;
        }
        
        if (loaded)
            freeData();

        loadData_tables(fi, table_count);
    }
    fclose(fi);

    updateMap();

    loaded=true;
}

void DataBase::updateMap(){
    name2tid.clear();
    for (size_t i=0;i<tables.size();i++)
        name2tid[tables[i].name]=i;
}

void DataBase::createNew(string _path, string _name){
    if (_name.size()>30){
        printf_error("error: database name is too long!\n");
        return;
    }
    if (loaded)
        freeData();
    
    path=_path;
    name=_name;
}

void DataBase::saveData(){
    string fullpath=path+name+".db";
    FILE *fo=fopen(fullpath.c_str(), "wb");
    if (!fo){
        printf_error("error: failed to save db");
        return;
    }

    char predata[64];
    strcpy(predata, "simplesql");
    strcpy(predata+10, "0.1");
    predata[63]=(char)((int)tables.size()-1)/2;
    fwrite(predata, 64, 1, fo);

    for (Table &tb: tables)
        tb.saveData(fo);
}
