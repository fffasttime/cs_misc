#include "storage.h"
#include <string.h>
#include <stdlib.h>

void Table::loadData(FILE *fi, string _name, int _item_count){
    assert(!loaded);
    name = _name;
    count = _item_count;

    data.resize(count);
    printf_debug("debug: loading table '%s', %d items\n", 
        name.c_str(), count);

    try{
        for (int i=0;i<count;i++){
            data[i] = (PRecord_t)malloc(field.length);
            if (data[i]==nullptr){
                printf_error("memory error: when creating table '%s', fail to allocate %d space", 
                    field.length);
                throw IOException();
            }
            if (!fread(data[i], field.length, 1, fi)){
                printf_error("error: when loading table '%s', unable to read file.\n", name.c_str());
                throw IOException();
            }
        }
    }
    catch(const IOException& e){
        loaded = true;
        freeData();
        throw IOException();
    }
    
    loaded=true;
}

void Table::freeData(){
    assert(loaded);
    for (const auto &p : data)
        if (p != nullptr)
            free(p);
    loaded=false;
}

void Table::saveData(FILE *fo){
    printf_debug("debug: saving table '%s'\n", name.c_str());
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
    if (fread(predata, 64, 1, fi) == 0){
        printf_error("error: bad db file format, unknown head\n");
        return false;
    }
    // check head
    if (strncmp(predata,"simplesql",9)){
        printf_error("error: bad db file format, maybe another dbms?\n");
        return false;
    }
    // check version
    const char cur_version[]="0.1";
    char file_version[10];
    strncpy(file_version, predata+10, 10);
    if (strncmp(file_version, "0.1", 9)){
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
    fields.emplace_back(FieldType::nchar,31,"name");
    fields.emplace_back(FieldType::int32,"count");
    fields.emplace_back(FieldType::int32,"count_field");
    return FieldInfo(fields);
}
FieldInfo getTableMetaFieldInfo(){
    vector<FieldCellInfo> fields;
    fields.emplace_back(FieldType::int32,"type");
    fields.emplace_back(FieldType::int32,"extra");
    fields.emplace_back(FieldType::nchar,31,"name");
    return FieldInfo(fields);
}

DataBase::~DataBase(){
    if (loaded){
        freeData();
    }
}

void DataBase::freeData(){
    assert(loaded);
    loaded=false;
    tables.clear();
    tbmeta.clear();
    dbmeta.freeData();
    name_tab.clear();
}

void DataBase::loadData_tables(FILE *fi, int table_count){

    //load dbmeta
    dbmeta.loadData(fi, "__dbmeta", table_count);

    // load table metadata
    for (int i=0;i<table_count;i++){
        tbmeta.emplace_back(getTableMetaFieldInfo());
        tbmeta.back().loadData(fi, string("__tbmeta_")+dbmeta.read(i,"name").nchar(), 
                        dbmeta.read(i,"count_field").int32());
    }

    // load ordinary data
    for (int i=0;i<table_count;i++){
        vector<FieldCellInfo> fields;
        
        Table &cur=tbmeta[i];
        for (int j=0;j<cur.count;j++)
            fields.emplace_back((FieldType)cur.read(j,"type").int32(),
                            cur.read(j,"extra").int32(), 
                            cur.read(j,"name").nchar());
        tables.emplace_back(FieldInfo(fields));
        tables.back().loadData(fi, dbmeta.read(i,"name").nchar(), 
                        dbmeta.read(i,"count").int32());
    }
}

/**
 * 'use <schema>' will load a exist database
 */
void DataBase::loadData(string _path, string _name){
    path=_path;
    name=_name;
    
    string fullpath=path+name+".db";
    
    FILE *fi=fopen(fullpath.c_str(), "rb");
    
    //check
    if (!fi){
        printf_error("error: can't open file %s\n",path.c_str());
        return;
    }
    int table_count;
    if (!checkDBFile(fi, table_count)){
        fclose(fi);
        return;
    }

    // start loading
    if (loaded) freeData();
    printf_debug("debug: schema '%s' load, %d tables\n", name.c_str(), table_count);

    try{
        loadData_tables(fi, table_count);
    }
    catch(const IOException& e){
        loaded = true; // for freeData()
        freeData();
        fclose(fi);
        return;
    }

    updateMap();
    loaded=true;
    
    fclose(fi);
}

void DataBase::updateMap(){
    name_tab.clear();
    name_tab[dbmeta.name] = &dbmeta;
    for (auto &tab : tables)
        name_tab[tab.name] = &tab;
    for (auto &tab : tbmeta)
        name_tab[tab.name] = &tab;
}

void DataBase::createNew(string _path, string _name){
    if (_name.size()>30){
        printf_error("error: database name is too long!\n");
        return;
    }
    if (loaded)
        freeData();
    
    dbmeta.name="__dbmeta";
    dbmeta.loaded = true; // don'n need to load file

    name_tab[dbmeta.name] = &dbmeta;
    loaded = true;
    path = _path;
    name = _name;
}

void DataBase::saveData(){
    string fullpath=path+name+".db";
    FILE *fo=fopen(fullpath.c_str(), "wb");
    if (!fo){
        printf_error("error: failed to save db");
        return;
    }
    
    printf_debug("debug: saving data of schema '%s'\n", name.c_str());
    char predata[64];
    strcpy(predata, "simplesql");
    strcpy(predata+10, "0.1");
    predata[63]=(char)tables.size();
    fwrite(predata, 64, 1, fo);

    dbmeta.saveData(fo);
    for (Table &tb: tbmeta)
        tb.saveData(fo);
    for (Table &tb: tables)
        tb.saveData(fo);
}
