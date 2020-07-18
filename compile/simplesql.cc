#include "simplesql.h"
#include <dirent.h> //file system
#include "storage.h"

DataBase db;

void hintCMD(){
    if (db.loaded)
        printf((db.name+"> ").c_str());
    else
    	printf("SQL> ");
}

vector<string> listDir(string path){
    std::vector<string> list;
    DIR *pdir;
    dirent *pdirent;
    if (!(d=opendir(path.c_str())))
        fatal("fatal: fail to open db dir");
    while ((pdirent=readdir(pdir))!=NULL){
        if (pdirent->d_type==8){
            list.emplace_back(pdirent->d_name);
        }
    }
    closedir(pdir);
}

/**
 * print inital information
 */
void initDB(){
    const string path("data/");
    auto files=list_dir(path);
    for (const auto &s:files){
        if (s.size()>3 && s.substr(s.size()-3)==".db"){
            printf_info("info: found db file %s\n",s.c_str());
        }
    }
    printf_info("type use <schema name> to load a db file");
}

void createTable(char *tableval, create_items_def *crtitem_root){
    printf("trying create %s\n", tableval);
}

void selection(item_def *item_first, table_def *table_first, conditions_def *con_root){
    
}
