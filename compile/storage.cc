#include "storage.h"
#include <string.h>
#include <stdlib.h>
#include <exception>

void Table::loadData(FILE *fi){
    assert(!loaded);
    data.resize(item_count);

    for (int i=0;i<item_count;i++){
        data[i]=malloc(field.length);
        fread(data[i],field.length,1,fi);
    }

    loaded=1;
}
void Table::freeData(){
    assert(loaded);
    for (const auto &p:data)
        if (p!=nullptr) 
            free(p);
}

Table::~Table(){
    for (const auto &p:data)
        if (p!=nullptr) 
            free(p);
}

/**
 * read first 64 bytes
 */
bool checkDBFile(File *fi, int &table_count){
    char* predata[64];
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
    vector<FieldCellInfo> fields;=
    fields.emplace_back(FieldType::nchar,32,"name");
    fields.emplace_back(FieldType::int32,"count");
    return FieldInfo(fields);
}
FieldInfo getTableMetaFieldInfo(){
    vector<FieldCellInfo> fields;
    fields.emplace_back(FieldType::int32,"type");
    fields.emplace_back(FieldType::int32,"extra");
    fields.emplace_back(FieldType::nchar,32,"name");
    return FieldInfo(fields);
}

void DataBase::loadData(string _path, string _name){
    path=_path;
    name=_name;
    
    string fullpath=path+name+".db";
    
    FILE *fi=fopen(fullpath.c_str(), "rb");
    if (!fi){
        printf_info("error: can't open file %s\n",path.c_str());
        return false;
    }

    int table_count;
    if (!check_dbfile(fi, table_count)){
        fclose(fi);
        return;
    }
    
    if (!loaded)
        loaded=true;

    tables.emplace_back(getDBMetaFieldInfo);
    tables[0].loadData(fi, "__dbmeta", table_count);

    for (int i=0;i<table_count;i++){
        auto fieldinfo=getTableMetaFieldInfo();
        tables.emplace_back(fieldinfo);
        
    }

    fclose(fi);
}
