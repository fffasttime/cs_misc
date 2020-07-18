/**
 * This header file defines some data structue used for sql parser
*/

#ifndef SIMPLESQL_H
#define SIMPLESQL_H

#include <stdio.h>
#include <string.h>
#include "util.h"

void db_init();

/*'CREATE TABLE table_id (create_items_def)'*/
struct create_items_def{ 
	/*field name*/
	char *field; 
	FieldType type;    
	create_items_def *next;
};

struct item_def 	/*item info when execute 'SELECT items_list FROM ...' etc.*/
{
	char *field;	/*name*/
	struct field *pos;	/*real pos*/
	struct item_def *next;
};

struct conditions_def /*codition binary tree node type*/
{
	/*INT:0  STRING:1*/	
	int type;  
	/*item*/
	item_def *litem; 
	int intv;		
	char *strv;
	/* '=':1 | '>':2 | '<':3 | '>=':4 | '<=':5 | '!=':6 | 'AND':7 | 'OR':8 */		
	int cmp_op;
	conditions_def *left;
	conditions_def *right;
};

struct table_def	/*'SELECT * FROM tabel_list WHERE ...'*/
{
	/*table_name*/
	char *table;
	/*real_table_pos*/
	Table *pos;
	table_def *next;
};

void createTable(char *tableval, create_items_def *crtitem_root);
void selection(item_def *item_first, table_def *table_first, conditions_def *con_root);

#endif
