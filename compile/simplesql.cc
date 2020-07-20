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

void createTable(char *name, create_item_def *crtitem){
    if (strlen(name)>30){
        printf_error("error: table name is too long!\n");
        return;
    }
    printf_debug("debug: trying create %s\n", name);

    //create fieldinfo
    vector<FieldCellInfo> fields;
    for (const auto &it : *crtitem){
        //printf_debug("  field %u, type=%d, ex=%d, name=%s\n", \
            fields.size(), it.type, it.extra, it.name);
        fields.emplace_back(it.type, it.extra, it.name);
        const auto &cur=fields.back();
        if (cur.name.size()>30){
            printf_error("error: one field name is too long!\n");
            return;
        }
        if (cur.type==FieldType::nchar){
            if (cur.extra<1){
                printf_error("error: CHAR(n) field length should larger than 0\n");
                return;
            }
            if (cur.extra>255){
                printf_error("error: CHAR(n) field length too large\n");
            }
        }
    }
    //insert into dbmeta
    vector<value_def_unit> values;
    values.push_back({name,FieldType::nchar}); //name
    values.push_back({0,FieldType::int32}); //count
    values.push_back({fields.size(),FieldType::int32}); //count_field
    insertRecord("__dbmeta",values,nullptr);
    
    //insert tbmeta
    db.tables.emplace_back(getTableMetaFieldInfo());
    db.tables.back().name=string("__tbmeta_")+name;
    for (const auto &it: fields){
        values.clear();
        
        insertRecord(cur.name.c_str(), values, nullptr);
    }

    //insert table
    tb.tables.emplace_back(FieldInfo(fields));
    db.tables.back().name=name;


}

void selection(select_item_def *item, 
                table_def *table, conditions_def *con_root){
    
}

/**
 * check weather a value matchs its field when insert or update
 * TODO: UNIQUE/PRIMARKEY/... check
 */
bool checkValueField(const value_def_unit &v, const FieldCellInfo &f){
    if (v.type!=f.type){
        printf_error("error: type mismatch on '%s'\n", f.name.c_str());
        return false;
    }
    if (v.type==FieldType::nchar){
        int len=strlen(v.value.strval);
        if (len>f.extra){
            printf_error("error: field '%s' is CHAR(%d), but input string has %d characters\n",
                f.name.c_str(), f.extra, len);
        }
    }
    return true;
}

void writeValueField(PRecord_t p, int pos, const Table& tb, const value_def_unit &v){
    auto p1=(void *)(p+tb.field.offset[pos]);
    switch (tb.field.fields[pos].type)
    {
    case FieldType::int32:
        *(int *)p1 = v.value.intval;
        break;
    case FieldType::nchar:
        strcpy((char *)p1, v.value.strval);
        break;
    }
}

/**
 * insert single record
 */
void insertRecord(const char *name, value_def *val, select_item_def *selitem, bool nocheck = false){
    if (!db.name2tid.count(name)){
        printf_error("error: table '%s' does not exist\n", name);
        return;
    }
    Table &tb=db.tables[db.name2tid[name]];
    if (nocheck) goto next;
    // check fields
    auto &fields = tb.field.fields;
    if (selitem==nullptr){
        if (val->size()!=tb.field.fields.size()){
            printf_error("error: table '%s' has %u fields, but %u in input value\n",
                    name, fields.size(), val->size());
            return;
        }
        for (size_t i=0;i<val->size();i++){
            if (!checkValueField((*val)[i], fields[i]))
                return;
        }
    }
    else{
        if (selitem->size()!=val->size()){
            printf("error: number mismatch, %u values but %u fields\n",
                val->size(), selitem->size());
        }
        for (size_t i=0;i<selitem->size();i++){
            auto it=tb.field.name2fid.find((*selitem)[i]);
            if (it==tb.field.name2fid.end()){
                printf_error("error: table '%s' has no field named '%s'\n", name, it);
                return;
            }
            if (!checkValueField(fields[it->second]))
                return;
        }
    }
    next:
    //generate record
    auto record=(PRecord_t)malloc(tb.field.length);
    if (selitem==nullptr){
        for (size_t i=0;i<val->size();i++)
            writeValueField(record, i, tb, (*val)[i]);
    }
    else{
        for (size_t i=0;i<selitem->size();i++)
            writeValueField(record, tb.field.name2fid[(*selitem)[i]], tb, (*val)[i]);
    }
    tb.data.push_back(record);
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
