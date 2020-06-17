#ifndef SIMPLESQL_H
#define SIMPLESQL_H

#include <stdio.h>
#include <string.h>

struct create_item_def{
	char *field;
	int type;
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
	item_def *litem; /*item*/
	int intv;		
	char *strv;		
	int cmp_op;		/* '=':1 | '>':2 | '<':3 | '>=':4 | '<=':5 | '!=':6 | 'AND':7 | 'OR':8 */
	conditions_def *left;
	conditions_def *right;
}

void createTable(char *tableval, create_item_def *crtitem_root);
void selection(item_def *item_first, table_def *table_first, conditions_def *con_root);

#endif
