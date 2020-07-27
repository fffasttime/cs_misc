/**
 * This header file defines operations and data structue used for sql parser
*/

#ifndef SIMPLESQL_H
#define SIMPLESQL_H

#include <stdio.h>
#include <string.h>
#include <vector>
#include <array>
#include "common.h"
#include "storage.h"
using std::vector;
class CommandException:exception{};

/* field info */
struct create_item_def_unit{ 
	FieldType type;
	int extra;  
	char *name;
	char constraint;
};
/* CREATE TABLE table_id (create_items_def) */
typedef vector<create_item_def_unit> create_item_def;

struct select_item_def_unit{ 
	char *name;
	char *tabname;
	char *rename;
};
/* SELECT select_item_def FROM ... */
typedef vector<select_item_def_unit> select_item_def;

union KeyValue{
	int intval;
	float floatval;
	const char *strval;
};

/* INSERT INTO table_name VALUES (value_def) */
struct value_def_unit{
	KeyValue value;
	FieldType type;
};
typedef vector<value_def_unit> value_def;

class Searcher;
/*codition binary tree node type*/
struct conditions_def{
	/* binary_op:0  NUMBER:1 STRING:2 ID:3 FLOAT:4 SELECTION:5*/	
	int type;  
	/** 
	 * NUMBER or cmp_op_id 
	 * '=':1 | '>':2 | '<':3 | '>=':4 | '<=':5 | '!=':6 | 'AND':7 | 'OR':8 | 'IN':9
	 */	
	int intv;
	/* STRING or ID*/
	char *strv;	
	double floatv;
	/* extra id */
	char *tablev;
	conditions_def *left;
	conditions_def *right;
	Searcher *subselect;
	string to_str();
};

/*SELECT * FROM tabel_list WHERE ...*/
typedef vector<string> table_def;

extern int output_mode;
void hintCMD();
vector<string> listDir(const char *path);
void initDB();

void createTable(char *name, create_item_def *crtitem);
void dropTable(char *name);

void showStmt(char *cmd, char *from);
void createDatabase(char *name);
void useDatabase(char *name);
void saveDatabase();
int insertRecord(const char *name, value_def *val, select_item_def *selitem, bool nocheck = false);
void insertRecord_user(const char *name, value_def *val, select_item_def *selitem);

void selection(select_item_def *item, table_def *table, conditions_def *con_root);
void deleteItem(char *name, conditions_def *con_root);
void updateItem(char *name, char *fname, value_def_unit val, conditions_def *con_root);

void exitSql();

/**
 * means the item position in table
 * -1 stands for any.
 */
typedef std::array<short, 6> ItemTuple;
typedef vector<ItemTuple> ItemSet;

class Searcher{
private:
	vector<Table *> tabs;
	select_item_def *projname;
	conditions_def *con_root;
	void debug(conditions_def *cur, int dep = 0);
	pair<int, int> findFieldName(char *name, char *tabname);
public:
	bool succeed;

	ItemSet result;
	
	Searcher(select_item_def *item, table_def *table, conditions_def *con_root);
	ItemSet search(conditions_def *cur);
	ItemSet ItemSetFill(ItemTuple used, ItemSet input);
	ItemSet conTableIn(conditions_def *left, conditions_def *right);
	ItemSet conLogic(conditions_def *cur);
	ItemSet conCompare(conditions_def *cur);
	ItemSet CompareTable0(conditions_def *left, conditions_def *right, int cmp_op);
	//left should be an ID
	ItemSet CompareTable1(conditions_def *left, conditions_def *right, int cmp_op);
	ItemSet CompareTable2(conditions_def *left, conditions_def *right, int cmp_op);

	void showResult();
	Table *doProjection();
};

#endif
