#include "simplesql.h"
#include <dirent.h> //file system
#include <iostream>
#include <algorithm>
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

value_def_unit getValDefUnit(int val, const FieldType &t){
    return (value_def_unit){{.intval=val}, t};
}
value_def_unit getValDefUnit(const char *val, const FieldType &t){
    return (value_def_unit){{.strval=val}, t};
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
        /* printf_debug("  field %u, type=%d, ex=%d, name=%s\n", 
            fields.size(), it.type, it.extra, it.name); */
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
    vector<value_def_unit> values = {
        getValDefUnit(name, FieldType::nchar),         //name
        getValDefUnit(0, FieldType::int32),            //count
        getValDefUnit((int)fields.size(), FieldType::int32) //count_field
    };
    insertRecord("__dbmeta", &values, nullptr);
    
    //insert tbmeta
    db.tbmeta.emplace_back(getTableMetaFieldInfo());
    db.tbmeta.back().name=string("__tbmeta_")+name;
    db.tbmeta.back().loaded=true;
    db.name_tab[db.tbmeta.back().name] = &db.tbmeta.back();
    for (const auto &it: fields){
        values = {
            getValDefUnit((int)it.type, FieldType::int32),        //name
            getValDefUnit(it.extra, FieldType::int32),       //count
            getValDefUnit(it.name.c_str(), FieldType::nchar) //count_field
        };
        insertRecord(db.tbmeta.back().name.c_str(), &values, nullptr);
    }

    //insert table
    db.tables.emplace_back(FieldInfo(fields));
    db.tables.back().name=name;
    db.tables.back().loaded=true;
    db.name_tab[name]=&db.tables.back();
}

/**
 * Check weather a value matchs its field when insert or update
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
void insertRecord(const char *name, value_def *val, select_item_def *selitem, bool nocheck){
    if (!db.name_tab.count(name)){
        printf_error("error: table '%s' does not exist\n", name);
        return;
    }
    Table &tb=*db.name_tab[name];
    auto &fields = tb.field.fields;
    if (nocheck) goto next;
    // check fields
    if (selitem==nullptr){
        if (val->size()!=tb.field.fields.size()){
            printf_error("error: table '%s' has %zu fields, but %zu in input value\n",
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
            printf("error: number mismatch, %zu values but %zu fields\n",
                val->size(), selitem->size());
        }
        for (size_t i=0;i<selitem->size();i++){
            auto it=tb.field.name2fid.find((*selitem)[i]);
            if (it==tb.field.name2fid.end()){
                printf_error("error: table '%s' has no field named '%s'\n", name, it);
                return;
            }
            if (!checkValueField((*val)[i], fields[it->second]))
                return;
        }
    }
    next:
    //generate record
    auto record=(PRecord_t)calloc(tb.field.length, 1);
    if (record==nullptr){
        printf_error("memory error: fail to allocate %d space", tb.field.length);
        return;
    }
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

string conditions_def::to_str(){
    if (type==0){
        switch (intv)
        {
        case 1: return "=";
        case 2: return ">";
        case 3: return "<";
        case 4: return ">=";
        case 5: return "<=";
        case 6: return "!=";
        case 7: return "AND";
        case 8: return "OR";
        }
    }
    else if (type == 1){
        return std::to_string(intv);
    }
    else{
        return strv;
    }
}

void Searcher::debug(conditions_def *cur, int dep){
    if (cur->left)
        debug(cur->left, dep + 1);
    for (int i = 0; i < dep * 2; i++)
        putchar(' ');
    puts(cur->to_str().c_str());
    if (cur->right)
        debug(cur->right, dep + 1);
}

Searcher::Searcher(select_item_def *item, table_def *table, conditions_def *con_root)
    :projname(item),con_root(con_root)                
{
    // list tables
    for (const auto &tabname: *table){
        auto it=db.name_tab.find(tabname);
        if (it==db.name_tab.end()){
            printf("error: no table named '%s'\n", tabname.c_str());
            return;
        }
        tabs.push_back(it->second);
    }
    debug(con_root);
    try{
        // start search
        result=conLogic(con_root);
    }
    catch(const CommandException& e){
        ;
    }
    if (projname==nullptr){
        //show
        for (const auto &it:result){
            putchar('|');
            for (size_t i=0;i<tabs.size();i++){
                auto tab=tabs[i];
                for (size_t j=0;j<tab->field.fields.size();j++){
                    int offset=tab->field.offset[j];
                    auto field=tab->field.fields[j];
                    if (field.type==FieldType::int32){
                        printf("%d", tab->readof(it[i], offset).int32());
                    }
                    else{
                        puts( tab->readof(it[i], offset).nchar());
                    }
                }
            }
        }
    }
}

bool compareint(int x, int y, int cmp_op){
    switch (cmp_op)
    {
    case 1: return x==y;
    case 2: return x>y;
    case 3: return x<y;
    case 4: return x>=y;
    case 5: return x<=y;
    case 6: return x!=y;
    }
}
bool comparestr(char *x, char *y, int cmp_op){
    int ret=strcmp(x, y);
    switch (cmp_op)
    {
    case 1: return ret==0;
    case 2: return ret>0;
    case 3: return ret<0;
    case 4: return ret>=0;
    case 5: return ret<=0;
    case 6: return ret;
    }
}

const ItemTuple ItemTuple_True={-1, -1, -1, -1, -1, -1};
const ItemSet ItemSet_True=ItemSet({ItemTuple_True});
const ItemSet ItemSet_False=ItemSet();

ItemSet Searcher::CompareTable0(conditions_def *left, conditions_def *right, int cmp_op){
    if (left->type==1 || right->type==1){
        int lx=left->type==1?left->intv:atoi(left->strv);
        int rx=right->type==1?right->intv:atoi(right->strv);
        return compareint(lx, rx, cmp_op)?ItemSet_True:ItemSet_False;
    }
    else{
        return comparestr(left->strv, right->strv, cmp_op)?ItemSet_True:ItemSet_False;
    }
}

/**
 * @return pair<int,int> : <table id in Searcher.tab, field id in real table>
 */
pair<int, int> Searcher::findFieldName(char *name){
    for (size_t i=0;i<tabs.size();i++){
        if (tabs[i]->field.name2fid.count(name))
            return {(int)i, tabs[i]->field.name2fid[name]};
    }
    printf_error("error: field name '%s' did not found\n", name);
    throw CommandException();
}

ItemSet Searcher::CompareTable1(conditions_def *left, conditions_def *right, int cmp_op){
    auto id=findFieldName(left->strv);
    int tid=id.first;
    Table &tab=*tabs[tid];
    auto &field=tab.field.fields[id.second];
    int offset=tab.field.offset[id.second];
    ItemSet result;
    // TODO : optimize, add type check
    //enum
    for (size_t i=0;i<tab.data.size();i++){
        if(field.type==FieldType::int32){
            int x=tab.readof(i, offset).int32();
            if (compareint(x,right->intv,cmp_op)){
                auto ins=ItemTuple_True;
                ins[tid]=i;
                result.emplace_back(std::move(ins));
            }
        }
        else{
            char *x=tab.readof(i, offset).nchar();
            if (comparestr(x,right->strv,cmp_op)){
                auto ins=ItemTuple_True;
                ins[tid]=i;
                result.emplace_back(std::move(ins));
            }
        }
    }
    return result;
}

ItemSet Searcher::CompareTable2(conditions_def *left, conditions_def *right, int cmp_op){
    auto id1=findFieldName(left->strv);
    int tid1=id1.first;
    Table &tab1=*tabs[tid1];
    auto &field1=tab1.field.fields[id1.second];
    int offset1=tab1.field.offset[id1.second];

    auto id2=findFieldName(right->strv);
    int tid2=id2.first;
    Table &tab2=*tabs[tid2];
    auto &field2=tab2.field.fields[id2.second];
    int offset2=tab2.field.offset[id2.second];

    ItemSet result;
    // TODO : optimize, add type check
    //enum
    for (size_t i=0;i<tab1.data.size();i++){
        if(field1.type==FieldType::int32){
            int x=tab1.readof(i, offset1).int32();
            for (int j=0;j<tab2.data.size();j++){
                int y=tab2.readof(i, offset2).int32();
                if (compareint(x,y,cmp_op)){
                    auto ins=ItemTuple_True;
                    ins[tid1]=i;
                    ins[tid2]=j;
                    result.emplace_back(std::move(ins));
                }
            }
        }
        else{
            char *x=tab1.readof(i, offset1).nchar();
            for (int j=0;j<tab2.data.size();j++){
                char *y=tab2.readof(i, offset2).nchar();
                if (comparestr(x,y,cmp_op)){
                    auto ins=ItemTuple_True;
                    ins[tid1]=i;
                    ins[tid2]=j;
                    result.emplace_back(std::move(ins));
                }
            }
        }
    }
    return result;
}

ItemSet Searcher::conCompare(conditions_def *cur){
    if (cur->type){
        printf("error: '%s' can't be a condition\n", cur->to_str().c_str());
        throw CommandException();
    }
    if (cur->left->type == 3 && cur->right->type == 3){
        return CompareTable2(cur->left, cur->right, cur->intv);
    }
    else if (cur->left->type !=3 && cur->right->type != 3){
        return CompareTable0(cur->left, cur->right, cur->intv);
    }
    else{
        if (cur->left->type == 3)
            return CompareTable1(cur->left, cur->right, cur->intv);
        else
            return CompareTable1(cur->right, cur->left, cur->intv);
    }
}

ItemSet Searcher::ItemSetFill(ItemTuple used, ItemSet input){
    if (input.size()==0) return input;
    ItemSet result;
    for (size_t i=0;i<tabs.size();i++){
        if (used[i]==-1 && input[0][i]!=-1){
            for (auto it : input){
                for (int j=0;j<tabs[i]->data.size();j++){
                    it[i]=j;
                    result.push_back(it);
                }
            }
        }
        swap(input, result);
        result.clear();
    }
    return input;
}

ItemSet Searcher::conLogic(conditions_def *cur){
    if (cur->type){
        printf("error: '%s' can't be a condition\n", cur->to_str().c_str());
        throw CommandException();
    }
    if (cur->intv==7 || cur->intv==8){ // AND/OR
        ItemSet l=conLogic(cur->left), r=conLogic(cur->right);
        l=ItemSetFill(ItemTuple_True, l);
        r=ItemSetFill(ItemTuple_True, r);
        // this method is easy, but very slow
        l.insert(l.end(), r.begin(), r.end());
        std::sort(l.begin(),l.end());
        auto it=std::unique(l.begin(),l.end());
        if (cur->intv==7) //intersection
            return ItemSet(l.begin(), it);
        else //union
            return ItemSet(it, l.end());
    }
    return conCompare(cur); //not logic binary
}

void selection(select_item_def *item, table_def *table, conditions_def *con_root){
    Searcher Searcher(item, table, con_root);
}
