/**
 * This header file defines operations and data structue used for sql parser
*/

#ifndef SIMPLESQL_H
#define SIMPLESQL_H

#include <stdio.h>
#include <string.h>
#include <vector>
#include "common.h"
#include "storage.h"
using std::vector;

struct create_item_def_unit{ 
	/*field name*/
	FieldType type;
	int extra;  
	char *name;
};
/*CREATE TABLE table_id (create_items_def)*/
typedef vector<create_item_def_unit> create_item_def;

/* SELECT select_item_def FROM ...*/
typedef vector<string> select_item_def;

union KeyValue{
	int intval;
	char *strval;
};

/* INSERT INTO table_name VALUES (value_def) */
struct value_def_unit{
	KeyValue value;
	FieldType type;
};
typedef vector<value_def_unit> value_def;

/*codition binary tree node type*/
struct conditions_def{
	/*INT:0  STRING:1*/	
	int type;  
	/*item*/
	// !-- TODO: check again
	char *litem; 
	int intv;		
	char *strv;
	/* '=':1 | '>':2 | '<':3 | '>=':4 | '<=':5 | '!=':6 | 'AND':7 | 'OR':8 */		
	int cmp_op;
	conditions_def *left;
	conditions_def *right;
};

/*SELECT * FROM tabel_list WHERE ...*/
typedef vector<string> table_def;

void hintCMD();
vector<string> listDir(const char *path);
void initDB();

void createTable(char *name, create_item_def *crtitem);
void selection(select_item_def *item, table_def *table, conditions_def *con_root);

void createDatabase(char *name);
void useDatabase(char *name);
void saveDatabase();
void insertRecord(char *name, value_def *val, select_item_def *selitem, bool nocheck);


#endif
