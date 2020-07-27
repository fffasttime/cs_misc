#include "simplesql.h"
#include <set>
#include <dirent.h> //file system
#include <iostream>
#include <algorithm>
DataBase db;

int output_mode;

void hintCMD(){
    if (output_mode&1) return;
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
    if (!db.loaded){
        printf_error("error: select a schema first");
        return;
    }

    if (strlen(name)>20){
        printf_error("error: table name is too long!\n");
        return;
    }

    //create fieldinfo
    vector<FieldCellInfo> fields;
    for (const auto &it : *crtitem){
        /* printf_debug("  field %u, type=%d, ex=%d, name=%s\n", 
            fields.size(), it.type, it.extra, it.name); */
        fields.emplace_back(it.type, it.extra, it.name, it.constraint);
        
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
    printf_debug("debug: creating table '%s' with %zu fields\n", name, fields.size());
    //insert into dbmeta
    vector<value_def_unit> values = {
        getValDefUnit(name, FieldType::nchar),         //name
        getValDefUnit(0, FieldType::int32),            //count
        getValDefUnit((int)fields.size(), FieldType::int32) //count_field
    };
    insertRecord("__dbmeta", &values, nullptr);
    
    // insert tbmeta
    db.tbmeta.emplace_back(getTableMetaFieldInfo());
    db.tbmeta.back().name=string("__tbmeta_")+name;
    db.tbmeta.back().loaded=true;
    db.updateMap();
    for (const auto &it: fields){
        values = {
            getValDefUnit((int)it.type, FieldType::int32), 
            getValDefUnit(it.extra, FieldType::int32),      
            getValDefUnit(it.name.c_str(), FieldType::nchar),
            getValDefUnit(it.constraint, FieldType::int8)
        };
        insertRecord(db.tbmeta.back().name.c_str(), &values, nullptr);
    }

    // insert table
    db.tables.emplace_back(FieldInfo(fields));
    db.tables.back().name=name;
    db.tables.back().loaded=true;
    db.updateMap();
}

void dropTable(char *name){
    if (db.name_tab.count(name)==0){
        printf_error("error: table '%s' does not exist\n", name);
        return;
    }
    printf_info("info: dropping table '%s'\n", name);
    for (size_t i=0;i<db.tables.size();i++)
        if(db.tables[i].name==name){
            // remove table
            db.tables[i].freeData();
            db.tables.erase(db.tables.begin()+i);
            db.tbmeta[i].freeData();
            db.tbmeta.erase(db.tbmeta.begin()+i);
            // delete meta item
            free(db.dbmeta.data[i]);
            db.dbmeta.data.erase(db.dbmeta.data.begin()+i);
            break;
        }
    db.updateMap();
}

void showStmt(char *cmd, char * from){
    if (string(cmd)=="databases"){
        initDB();
    }
    else if (string(cmd)=="columns"){
        if (from==nullptr){
            printf_error("syntax error, expecting FORM\n");
            return;
        }
        if (!db.name_tab.count(from)){
            printf_error("error: table '%s' does not exist\n", from);
            return;
        }
        table_def tdef({string("__tbmeta_")+from});
        selection(nullptr, &tdef, nullptr);
    }
    else if (string(cmd)=="tables"){
        table_def tdef({"__dbmeta"});
        selection(nullptr, &tdef, nullptr);
    }
    else{
        printf_error("syntax error, unexpected ID, expecting TABLES, DATABASES or COLUMNS\n");
    }
}

/**
 * Check weather a value matchs its field when insert or update
 * TODO: UNIQUE/PRIMARKEY/... check
 */
bool checkValueField(const value_def_unit &v, const FieldCellInfo &f){
    if (v.type==FieldType::Null && !f.allowNULL()){
        printf_error("error: field '%s' can't be NULL\n", f.name.c_str());
        return false;
    }
    if (f.type==FieldType::Float && v.type==FieldType::nchar){
        printf_error("error: type mismatch on '%s'\n", f.name.c_str());
        return false;
    }
    if (f.type==FieldType::nchar && v.type!=FieldType::nchar){
        printf_error("error: type mismatch on '%s'\n", f.name.c_str());
        return false;
    }
    if (v.type==FieldType::nchar){
        int len=strlen(v.value.strval);
        if (len>f.extra){
            printf_error("error: field '%s' is CHAR(%d), but input string has %d characters\n",
                f.name.c_str(), f.extra, len);
            return false;
        }
    }
    return true;
}

void writeValueField(PRecord_t p, int pos, const Table& tb, const value_def_unit &v){
    auto p1=(void *)(p+tb.field.offset[pos]);
    if (v.type==FieldType::Null){
        *((char *)p1+tb.field.fields[pos].length())=1;
        return;
    }
    switch (tb.field.fields[pos].type)
    {
    case FieldType::int32:
        *(int *)p1 = v.value.intval;
        break;
    case FieldType::int8:
        *(char *)p1 = v.value.intval;
        break;
    case FieldType::nchar:
        strcpy((char *)p1, v.value.strval);
        break;
    case FieldType::Float:{
        *(float *)p1 = (float)(v.type==FieldType::Float ? v.value.floatval:v.value.intval);
        break;
    }
    default:
        fatal("fatal: internal error in writeValueField(), type=%d\n",
            (int)tb.field.fields[pos].type);
    }
}

bool checkUnique(int pos, Table &tb, const value_def_unit &v){
    if ((tb.field.fields[pos].constraint&2)==0) // allow
        return true;
    int offset=tb.field.offset[pos];
    FieldCellInfo &f=tb.field.fields[pos];
    if (f.type==FieldType::int32){
        int y=v.value.intval;
        for (size_t i=0;i<tb.data.size();i++){
            if (tb.readof(i, offset).int32()==y){
                printf_error("error: field '%s' should be UNIQUE\n", f.name.c_str());
                return false;
            }
        }
    }
    else if (f.type==FieldType::nchar){
        for (size_t i=0;i<tb.data.size();i++)
            if (strcmp(tb.readof(i, offset).nchar(),v.value.strval)==0){
                printf_error("error: field '%s' should be UNIQUE\n", f.name.c_str());
                return false;
            }
    }
    return true;
}

/**
 * insert single record
 */
int insertRecord(const char *name, value_def *val, select_item_def *selitem, bool nocheck){
    if (!db.name_tab.count(name)){
        printf_error("error: table '%s' does not exist\n", name);
        return 1;
    }
    Table &tb=*db.name_tab[name];
    auto &fields = tb.field.fields;
    value_def_unit nullunit;
    nullunit.type=FieldType::Null;

    if (nocheck) goto next;
    // check fields
    if (selitem==nullptr){
        if (val->size()!=tb.field.fields.size()){
            printf_error("error: table '%s' has %zu fields, but %zu in input value\n",
                    name, fields.size(), val->size());
            return 1;
        }
        for (size_t i=0;i<val->size();i++){
            if (!checkValueField((*val)[i], fields[i]))
                return 1;
            if (!checkUnique(i, tb, (*val)[i]))
                return 1;
        }
    }
    else{
        if (selitem->size()!=val->size()){
            printf_error("error: number mismatch, %zu values but %zu fields\n",
                val->size(), selitem->size());
        }
        vector<bool> vis(tb.field.fields.size(), false);
        for (size_t i=0;i<selitem->size();i++){
            auto it=tb.field.name2fid.find((*selitem)[i]->name);
            if (it==tb.field.name2fid.end()){
                printf_error("error: table '%s' has no field named '%s'\n", name, (*selitem)[i]);
                return 1;
            }
            if (!checkValueField((*val)[i], fields[it->second]))
                return 1;
            if (!checkUnique(it->second, tb, (*val)[i]))
                return 1;
            vis[it->second]=true;
        }
        // check null
        for (size_t i=0;i<vis.size();i++)
            if(vis[i]==false)
                if (!checkValueField(nullunit,fields[i]))
                    return 1;
    }
    next:
    //generate record
    auto record=(PRecord_t)calloc(tb.field.length, 1);
    if (record==nullptr){
        printf_error("memory error: fail to allocate %d space", tb.field.length);
        return 1;
    }
    if (selitem==nullptr){
        for (size_t i=0;i<val->size();i++)
            writeValueField(record, i, tb, (*val)[i]);
    }
    else{
        vector<bool> vis(tb.field.fields.size(), false);
        for (size_t i=0;i<selitem->size();i++){
            int fid=tb.field.name2fid[(*selitem)[i]->name];
            writeValueField(record, fid, tb, (*val)[i]);
            vis[fid]=true;
        }
        // check null
        for (size_t i=0;i<vis.size();i++)
            if(vis[i]==false)
                writeValueField(record, i, tb, nullunit);
    }
    tb.data.push_back(record);
    return 0;
}

void insertRecord_user(const char *name, value_def *val, select_item_def *selitem){
    if (name[0]=='_'){
        printf_error("error: can't change metadata by insert\n");
        return;
    }
    if (insertRecord(name, val, selitem))
        return; //fail
    //update count
    for (size_t i=0;i<db.tables.size();i++)
        if (db.tables[i].name==name){
            db.dbmeta.read(i,"count").int32()++;
            break;
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

void exitSql(){
    if (db.loaded){
        printf_info("info: saving data\n");
        saveDatabase();
    }
    printf_info("info: all data saved, exiting\n");
    exit(0);
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
        case 9: return "IN";
        default: fatal("fatal: unknown op %d in function %s\n", intv, __func__);
        }
    }
    else if (type == 1){
        return std::to_string(intv);
    }
    else if (type == 2){
        return string("'")+strv+"'";
    }
    else if (type == 4){
        return std::to_string(floatv);
    }
    else if (type == 3){
        if (tablev!=nullptr)
            return string(tablev)+"."+strv;
        return strv;
    }
    else if (type == 5){
        return "<sub selection>";
    }
    return "error-type";
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
    succeed=0;
    // list tables
    for (const auto &tabname: *table){
        auto it=db.name_tab.find(tabname);
        if (it==db.name_tab.end()){
            printf_error("error: no table named '%s'\n", tabname.c_str());
            return;
        }
        tabs.push_back(it->second);
    }
    if (con_root!=nullptr)
        debug(con_root);
    try{
        // start search
        result=search(con_root);
    }
    catch(const CommandException& e){
        return;
    }
    succeed=1;
}

Table *Searcher::doProjection(){
    try{
        vector<FieldCellInfo> fields;
        vector<pair<int, int>> fids; 
        if (projname!=nullptr){
            fids.reserve(projname->size());
            for (auto &it : *projname){
                auto res=findFieldName(it->name, it->tabname);
                auto &tab=tabs[res.first];
                auto &fid=tab->field.fields[res.second];
                fids.emplace_back(res.first, res.second);
                if (it->rename==nullptr)
                    fields.push_back(FieldCellInfo(fid.type, fid.extra, fid.name));
                else
                    fields.push_back(FieldCellInfo(fid.type, fid.extra, it->rename));
            }
        }
        else{
            int i=0;
            for (auto &tab:tabs){
                for (size_t j=0;j<tab->field.fields.size();j++){
                    auto &fid=tab->field.fields[j];
                    fids.emplace_back(i,j);
                    fields.push_back(FieldCellInfo(fid.type, fid.extra, fid.name));
                }
                i++;
            }
        }
        Table *result_table=new Table(fields);
        auto &tab=*result_table;
        tab.loaded=true;
        tab.data.reserve(result.size());
        // copy data to new table
        // TODO : redesign it
        for (auto &it:result){
            auto line=(PRecord_t)malloc(tab.field.length);
            for (size_t i=0;i<fids.size();i++)
                memcpy(line+tab.field.offset[i], 
                tabs[fids[i].first]->data[it[fids[i].first]]+
                    tabs[fids[i].first]->field.offset[fids[i].second],
                tab.field.fields[i].length()+1);
            tab.data.push_back(line);
        }
        return result_table;
    }
    catch (CommandException &e){
        return nullptr;
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
    default: ;
    }
    fatal("fatal: interneral error, cmp_op=%d in function compareint\n", cmp_op);
}
bool comparefloat(float x, float y, int cmp_op){
    switch (cmp_op)
    {
    case 1: return x==y;
    case 2: return x>y;
    case 3: return x<y;
    case 4: return x>=y;
    case 5: return x<=y;
    case 6: return x!=y;
    default: ;
    }
    fatal("fatal: interneral error, cmp_op=%d in function comparefloat\n", cmp_op);
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
    fatal("fatal: interneral error, cmp_op=%d in function comparestr\n", cmp_op);
}

const ItemTuple ItemTuple_True={-1, -1, -1, -1, -1, -1};
const ItemSet ItemSet_True=ItemSet({ItemTuple_True});
const ItemSet ItemSet_False=ItemSet();

ItemSet Searcher::search(conditions_def *con_root){
    if (con_root==nullptr) // no where clause
        return ItemSetFill(ItemTuple_True, ItemSet({ItemTuple_True}));
    return ItemSetFill(ItemTuple_True, conLogic(con_root));
}

ItemSet Searcher::CompareTable0(conditions_def *left, conditions_def *right, int cmp_op){
    /** TODO:finish compare */
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
pair<int, int> Searcher::findFieldName(char *name, char *tabname){
    if (tabname!=nullptr){
        for (size_t i=0;i<tabs.size();i++)
            if (tabs[i]->name==tabname){
                if (tabs[i]->field.name2fid.count(name))
                    return {(int)i, tabs[i]->field.name2fid[name]};
                printf_error("error: in WHERE clause, table '%s' has no field named '%s'\n",
                    tabname, name);
                throw CommandException();
            }
        printf_error("error: in WHERE clause, table '%s' does not exist\n", tabname);
        throw CommandException();
    }
    int cnt=0;
    for (size_t i=0;i<tabs.size();i++)
        if (tabs[i]->field.name2fid.count(name))
            cnt++;
    if (cnt>1){
        printf_error("error: field name '%s' is conflict between tables\n", name);
        throw CommandException();
    }
    for (size_t i=0;i<tabs.size();i++){
        if (tabs[i]->field.name2fid.count(name))
            return {(int)i, tabs[i]->field.name2fid[name]};
    }
    printf_error("error: field name '%s' did not found\n", name);
    throw CommandException();
}

ItemSet Searcher::CompareTable1(conditions_def *left, conditions_def *right, int cmp_op){
    auto id=findFieldName(left->strv, left->tablev);
    int tid=id.first;
    Table &tab=*tabs[tid];
    auto &field=tab.field.fields[id.second];
    int offset=tab.field.offset[id.second];
    ItemSet result;
    // TODO : optimize, add type check
    //enum
    for (size_t i=0;i<tab.data.size();i++){
        if(field.type==FieldType::int32 || field.type==FieldType::Float){
            float x;
            if (field.type==FieldType::int32) x=tab.readof(i, offset).int32();
            else x=tab.readof(i, offset).float32();
            tab.readof(i, offset).int32();
            float y;
            if (right->type==1)
                y=right->intv;
            else if (right->type==4)
                y=right->floatv;
            else{
                printf_error("error: compare type mismatch\n");
                throw CommandException();
            }
            if (comparefloat(x, y,cmp_op)){
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
    auto id1=findFieldName(left->strv, left->tablev);
    int tid1=id1.first;
    Table &tab1=*tabs[tid1];
    auto &field1=tab1.field.fields[id1.second];
    int offset1=tab1.field.offset[id1.second];

    auto id2=findFieldName(right->strv, right->tablev);
    int tid2=id2.first;
    Table &tab2=*tabs[tid2];
    auto &field2=tab2.field.fields[id2.second];
    int offset2=tab2.field.offset[id2.second];

    ItemSet result;
    // TODO : optimize, add type check
    //enum
    for (size_t i=0;i<tab1.data.size();i++){
        if(field1.type==FieldType::int32 || field1.type==FieldType::Float){
            float x;
            if (field1.type==FieldType::int32) x=tab1.readof(i, offset1).int32();
            else x=tab1.readof(i, offset1).float32();

            for (size_t j=0;j<tab2.data.size();j++){
                float y;
                if (field2.type==FieldType::int32) y=tab2.readof(i, offset2).int32();
                else if (field2.type==FieldType::Float) y=tab2.readof(i, offset2).float32();
                else{
                    printf_error("error: compare type mismatch\n");
                    throw CommandException();
                }
                if (comparefloat(x,y,cmp_op)){
                    auto ins=ItemTuple_True;
                    ins[tid1]=i;
                    ins[tid2]=j;
                    result.emplace_back(std::move(ins));
                }
            }
        }
        else{
            if (field2.type!=FieldType::nchar){
                printf_error("error: compare type mismatch\n");
                throw CommandException();
            }
            char *x=tab1.readof(i, offset1).nchar();
            for (size_t j=0;j<tab2.data.size();j++){
                char *y=tab2.readof(j, offset2).nchar();
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

/**
 * @brief 'WHERE <field> IN (sub_selection)' process
*/
ItemSet Searcher::conTableIn(conditions_def *left, conditions_def *right){
    if (left->type!=3){
        printf_error("error: expecting a field name before IN\n");
        throw CommandException();
    }
    if (right->type!=5){
        printf_error("error: expecting a table after IN\n");
        throw CommandException();
    }
    auto searcher=right->subselect;
    if (!searcher->succeed){
        free(searcher);
        throw CommandException();
    }
    auto tab2=searcher->doProjection();
    if (tab2==nullptr){
        free(searcher);
        throw CommandException();
    }
    auto id=findFieldName(left->strv, left->tablev);
    int tid=id.first;
    Table &tab=*tabs[tid];
    auto &field=tab.field.fields[id.second];
    int offset=tab.field.offset[id.second];

    ItemSet result;
    try{
    
    if(tab2->field.fields.size()>1){
        printf_error("error: sub selection result has more than 1 column\n");
        throw CommandException();
    }
    if(tab2->field.fields[0].type!=field.type){
        printf_error("error: field type mismatch between '%s' IN <sub selection>\n",
            field.name.c_str());
        
        throw CommandException();
    }
    if(field.type==FieldType::nchar){
        std::set<string> cur;
        for (size_t i=0;i<tab2->data.size();i++)
            cur.emplace(tab2->readof(i,0).nchar());
        for (size_t i=0;i<tab.data.size();i++)
            if (cur.count(tab.readof(i,offset).nchar())){
                auto ins=ItemTuple_True;
                ins[tid]=i;
                result.emplace_back(ins);
            }
    }
    else if(field.type==FieldType::Float){
        std::set<float> cur;
        for (size_t i=0;i<tab2->data.size();i++)
            cur.insert(tab2->readof(i,0).float32());
        for (size_t i=0;i<tab.data.size();i++)
            if (cur.count(tab.readof(i,offset).float32())){
                auto ins=ItemTuple_True;
                ins[tid]=i;
                result.emplace_back(ins);
            }
    }
    else{
        std::set<int> cur;
        for (size_t i=0;i<tab2->data.size();i++)
            cur.insert(tab2->readof(i,0).int32());
        for (size_t i=0;i<tab.data.size();i++)
            if (cur.count(tab.readof(i,offset).int32())){
                auto ins=ItemTuple_True;
                ins[tid]=i;
                result.emplace_back(ins);
            }
    }

    }
    catch (CommandException &e){
        tab2->freeData();
        free(tab2);
        free(searcher);
        throw e;
    }

    tab2->freeData();
    free(tab2);
    free(searcher);
    return result;
}

ItemSet Searcher::conCompare(conditions_def *cur){
    if (cur->type){
        printf_error("error: '%s' can't be a condition\n", cur->to_str().c_str());
        throw CommandException();
    }
    if (cur->intv==9){
        return conTableIn(cur->left, cur->right);
    }
    if (cur->right->type==5 || cur->left->type==5){
        printf_error("syntax error: SELECT clause can't be here\n");
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
        if (used[i]==-1 && input[0][i]==-1){
            for (auto it : input){
                for (size_t j=0;j<tabs[i]->data.size();j++){
                    it[i]=j;
                    result.push_back(it);
                }
            }
            swap(input, result);
            result.clear();
        }
    }
    return input;
}

ItemSet Searcher::conLogic(conditions_def *cur){
    if (cur->type){
        printf_error("error: '%s' can't be a condition\n", cur->to_str().c_str());
        throw CommandException();
    }
    if (cur->intv==7 || cur->intv==8){ // AND/OR
        ItemSet l=conLogic(cur->left), r=conLogic(cur->right);
        l=ItemSetFill(ItemTuple_True, l);
        r=ItemSetFill(ItemTuple_True, r);
        std::sort(l.begin(),l.end());
        std::sort(r.begin(),r.end());
        ItemSet result;
        if (cur->intv==8) //union
            std::set_union(l.begin(),l.end(),r.begin(),r.end(), std::back_inserter(result));
        else //intersection
            std::set_intersection(l.begin(),l.end(),r.begin(),r.end(), std::back_inserter(result));
        auto it=std::unique(result.begin(),result.end());
        result.erase(it,result.end());
        return result;
    }
    return conCompare(cur); //not logic binary
}

void selection(select_item_def *item, table_def *table, conditions_def *con_root){
    Searcher searcher(item, table, con_root);
    if (searcher.succeed){
        Table *tab = searcher.doProjection();
        if (tab!=nullptr){
            tab->showData();
            tab->freeData();
            free(tab);
        }
    }
}

void deleteItem(char *name, conditions_def *con_root){
    if (!db.name_tab.count(name)){
        printf_error("error: table '%s' does not exist\n", name);
        return;
    }
    auto& tab=*db.name_tab[name];
    if (con_root==nullptr){ //remove all
        tab.freeData();
        tab.loaded=true;
        printf_info("info: deleting all %zu items\n", tab.data.size());
        tab.data.clear();
        //update count
        for (size_t i=0;i<db.tables.size();i++)
            if (db.tables[i].name==name){
                db.dbmeta.read(i,"count").int32()=0;
                break;
            }
        return;
    }
    vector<string> tabname({name});
    Searcher searcher(nullptr, &tabname, con_root);

    vector<PRecord_t> newdata;
    vector<bool> removed;
    removed.resize(tab.data.size());
    for (auto& it:searcher.result)
        removed[it[0]]=1;
    for (size_t i=0;i<tab.data.size();i++){
        if (removed[i])
            free(tab.data[i]);
        else
            newdata.push_back(tab.data[i]);
    }
    printf_info("info: deleted %zu items, remain %zu items\n",
        tab.data.size()-newdata.size(), newdata.size());

    //update count
    for (size_t i=0;i<db.tables.size();i++)
        if (db.tables[i].name==name){
            db.dbmeta.read(i,"count").int32()=newdata.size();
            break;
        }

    newdata.swap(tab.data);
}

void updateItem(char *name, char *fname, value_def_unit val, conditions_def *con_root){
    if (!db.name_tab.count(name)){
        printf_error("error: table '%s' does not exist\n", name);
        return;
    }
    auto& tab=*db.name_tab[name];
    if (!tab.field.name2fid.count(fname)){
        printf_error("error: field '%s' does not exist in table '%s'\n", fname, name);
        return;
    }
    vector<string> tabname({name});
    Searcher searcher(nullptr, &tabname, con_root);
    int fid=tab.field.name2fid[fname];
    int offset=tab.field.offset[fid];

    if (!checkValueField(val, tab.field.fields[fid]))
        return;

    printf_info("info: updating %zu items\n", searcher.result.size());

    if (tab.field.fields[fid].type==FieldType::int32){
        for (auto& it:searcher.result)
            tab.readof(it[0], offset).int32()=val.value.intval;
    }
    else{
        for (auto& it:searcher.result)
            strcpy(tab.readof(it[0],offset).nchar(), val.value.strval);
    }
}
