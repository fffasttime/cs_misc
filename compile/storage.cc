#include "storage.h"
#include <string.h>
#include <stdlib.h>


void Table::loaddata(FILE *fi){
    assert(!loaded);
    data.resize(item_count);

    for (int i=0;i<item_count;i++){
        data[i]=malloc(field.length);
        fread(data[i],field.length,1,fi);
    }

    loaded=1;
}
void Table::freedata(){
    assert(loaded);
    for (const auto &p:data)
        if (p!=nullptr) 
            free(p);
}

Table::Table(/* args */){
    loaded=false;
    
}

Table::~Table(){
    for (const auto &p:data)
        if (p!=nullptr) 
            free(p);
}

DataBase::Database(string path, string dbname){
    FILE *fi=fopen((path+dbname).c_str(),"rb");
    char head_data[64];
    fread(head_data,128,1,);
    if (strcmp(head_data+10,"0.1.0")){

    }
}
