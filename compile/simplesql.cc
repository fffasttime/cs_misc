#include "simplesql.h"
#include <dirent.h> //file system
#include <iostream>
DataBase db;

void hintCMD(){
    if (db.loaded)
        std::cout<<db.name+"> ";
    else
        std::cout<<"sql> ";
}

vector<string> listDir(const char *path){
    vector<string> list;
    DIR *pdir;
    dirent *pdirent;
    if (!(pdir=opendir(path)))
        fatal("fatal: fail to open db dir %s\n", path);
    printf_debug("debug: listing folder, pdir=%p\n", pdir);
    while ((pdirent=readdir(pdir))!=NULL){
        //list all file
        if (pdirent->d_type==8)
            list.emplace_back(pdirent->d_name);
    }
    closedir(pdir);
    return list;
}

/**
 * print inital information
 */
void initDB(){
    vector<string> files=listDir(DB_DATA_PATH);
    for (const auto &s:files){
        if (s.size()>3 && s.substr(s.size()-3)==".db"){
            printf_info("info: found db file '%s'\n",s.c_str());
        }
    }
    printf_info("info: type 'use <schema name>' or create a new db\n");
}

void createTable(char *name, create_item_def *crtitem_begin){
    if (strlen(name)>30){
        printf_error("error: table name is too long!\n");
        return;
    }
    printf_debug("debug: trying create %s\n", name);
    
}

void selection(select_item_def *item_begin, table_def *table_begin, conditions_def *con_root){
    
}

/**
 * insert single record
 */
void insertRecord(char *name, value_def val_begin, select_item_def *selitem_begin){
    if (!db.name2tid.count(name)){
        printf("error: table '%s' does not exist\n", name);
        return;
    }
    if (selitem_begin==nullptr){
        
    }
    else{
        
    }
}

void createDatabase(char *name){
    db.createNew(DB_DATA_PATH, name);
}

void useDatabase(char *name){
    db.loadData(DB_DATA_PATH, name);
}

void saveDatabase(){
    db.saveData();
}
