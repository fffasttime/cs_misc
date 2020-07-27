%{

#include "simplesql.h"

int yylex();
int yyparse();
void yyerror(const char *str){
	fprintf(stderr, "%s\n", str);
}
extern "C" int yywrap(){
	return 1;
}

int main(int argc, char **argv){
	output_mode=0;
	initDB();
	if (argc>1){
		freopen(argv[1],"r",stdin);
		output_mode|=1;
	}
	hintCMD();
	yyparse();
	return 0;
}

%}

%union{
	int intval;
	char *strval;
	double floatval;
	create_item_def *crtitemval;
	create_item_def_unit crtitemvalu;
	select_item_def *selitemval;
	select_item_def_unit selitenvalu;
	value_def *valdef;
	value_def_unit valdefu;
	conditions_def *conval;
	table_def *tbval;
	Searcher *subselectsqlp;
}

%error-verbose

%token AND AS CREATE CHAR DATABASE DELETE DROP FLOAT EXIT FROM INSERT IN INT INTO 
%token NOT NULLV SAVE SET SCHEMA SELECT SHOW TABLE UNIQUE UPDATE USE VALUES WHERE
%token <floatval> FLOATVAL
%token <intval> NUMBER
%token <strval> STRING ID

/* binding content */
%type <intval> comparator constraint
%type <crtitemvalu> create_item
%type <crtitemval> create_items
%type <conval> condition conditions condition_item
%type <selitemval> items items_or_star
%type <selitenvalu> item
%type <tbval> tables
%type <valdefu> value
%type <valdef> value_list
%type <subselectsqlp> subselectsql

%left OR
%left AND

%%

statements: statements statementline | statementline
statementline: statement ';' 
			| statement '\n' {hintCMD();} 
			| error '\n' { /*error recovery*/
				yyerrok;
				hintCMD();
			}

statement: createsql | selectsql | exitsql | insertsql | dropsql | deletesql | updatesql
	| usestmt | savestmt | showstmt | %empty 

usestmt: USE ID {
		useDatabase($2);
	}
savestmt: SAVE{
		saveDatabase();
	}

showstmt: SHOW ID{
		showStmt($2, nullptr);
	}
	| SHOW ID FROM ID{
		showStmt($2, $4);
	}

selectsql: SELECT items_or_star FROM tables {
		selection($2, $4, nullptr);
	}
	| SELECT items_or_star FROM tables WHERE conditions {
		selection($2, $4, $6);
	}

items_or_star: items { $$ = $1; }
	| '*' {$$ = nullptr; }

exitsql: EXIT {
		exitSql();
	}

createsql: CREATE TABLE ID '(' create_items ')' {
		createTable($3, $5);
		delete $5;
	}
	| CREATE DATABASE ID{
		createDatabase($3);
	}
	| CREATE SCHEMA ID{
		createDatabase($3);
	}

constraint: %empty{
		$$=0;
	}
	| NOT NULLV{
		$$=1;
	}
	| UNIQUE{
		$$=2;
	}

create_item: ID INT constraint{
		$$.name=$1;
		$$.type=FieldType::int32;
		$$.constraint=$3;
	}
	| ID CHAR '(' NUMBER ')' constraint{
		$$.name=$1;
		$$.type=FieldType::nchar;
		$$.extra=$4;
		$$.constraint=$6;
	}
	| ID FLOAT constraint{
		$$.name=$1;
		$$.type=FieldType::Float;
		$$.constraint=$3;
	}

create_items: create_item { 
		$$ = new create_item_def({$1});
	}
	| create_items ',' create_item{
		$$ = $1;
		$$->emplace_back($3);
	}

dropsql: DROP TABLE ID {
		dropTable($3);
	}

deletesql: DELETE FROM ID {
		deleteItem($3, nullptr);
	}
	| DELETE '*' FROM ID{
		deleteItem($4, nullptr);
	}
	| DELETE FROM ID WHERE conditions{
		deleteItem($3, $5);
	}

insertsql: INSERT INTO ID VALUES '(' value_list ')'{
		insertRecord_user($3, $6, nullptr);
		delete $6;
	}
	| INSERT INTO ID '(' items ')' VALUES '(' value_list ')'{
		insertRecord_user($3, $9, $5);
		delete $5;
		delete $9;
	}

updatesql: UPDATE ID SET ID '=' value{
		updateItem($2, $4, $6, nullptr);
	}
	| UPDATE ID SET ID '=' value WHERE conditions{
		updateItem($2, $4, $6, $8);
	}

value: NUMBER {
		$$.value.intval = $1;
		$$.type = FieldType::int32;
	}
	| STRING {
		$$.value.strval = $1;
		$$.type = FieldType::nchar;
	}
	| FLOATVAL {
		$$.value.floatval = $1;
		$$.type = FieldType::Float;
	}
	| NULLV {
		$$.type = FieldType::Null;
	}

value_list: value { 
		$$ = new value_def({$1});
	}
	| value_list ',' value {
		$$ = $1;
		$$->emplace_back($3);
	}

tables:	ID {
		$$ = new table_def({$1});
	}
	| tables ',' ID {
		$$ = $1;
		$$->emplace_back($3);			
	}

item: ID{
		$$.tabname = nullptr;
		$$.name = $1;
		$$.rename = nullptr;
	}
	| ID '.' ID{
		$$.tabname = $1;
		$$.name = $3;
		$$.rename = nullptr;
	}
	| ID '.' ID AS STRING{
		$$.tabname = $1;
		$$.name = $3;
		$$.rename = $5;
	}

items: item { $$ = new select_item_def({$1}); }
	| items ',' item { 
		$$ = $1;
		$$->emplace_back($3);
	}

comparator:		  
	  '='     {$$ = 1;}
	| '>'     {$$ = 2;}
	| '<'     {$$ = 3;}
	| '>' '=' {$$ = 4;}
	| '<' '=' {$$ = 5;}
	| '<' '>' {$$ = 6;}
	| '!' '=' {$$ = 6;}
	| IN      {$$ = 7;}

subselectsql: SELECT items_or_star FROM tables {
		$$ = new Searcher($2, $4, nullptr);
	}
	| SELECT items_or_star FROM tables WHERE conditions {
		$$ = new Searcher($2, $4, $6);
	}

condition_item: ID {
		$$ = new conditions_def;
		$$->type = 3;
		$$->strv = $1;
		$$->tablev = nullptr;
		$$->left = $$->right = nullptr;
	}
	| ID '.' ID{
		$$ = new conditions_def;
		$$->type = 3;
		$$->strv = $3;
		$$->tablev = $1;
		$$->left = $$->right = nullptr;
	}
	| NUMBER {
		$$ = new conditions_def;
		$$->type = 1;
		$$->intv = $1;
		$$->left = $$->right = nullptr;
	}
	| STRING {
		$$ = new conditions_def;
		$$->type = 2;
		$$->strv = $1;
		$$->left = $$->right = nullptr;
	}
	| FLOATVAL {
		$$ = new conditions_def;
		$$->type = 4;
		$$->floatv = $1;
		$$->left = $$->right = nullptr;
	}
	| '(' subselectsql ')' {
		$$->type = 5;
		$$->subselect = $2;
		$$->left = $$->right = nullptr;
	}

condition: 
	condition_item comparator condition_item {
		$$ = new conditions_def;
		$$->type = 0;
		$$->intv = $2;
		$$->left = $1;
		$$->right = $3;	
	}

conditions: condition { $$ = $1; }
	|   '(' conditions ')' { $$ = $2; }
	|   conditions OR conditions {
		$$ = new conditions_def;
		$$->type = 0;
		$$->intv = 8;
		$$->left = $1;
		$$->right = $3;	
	}
	|   conditions AND conditions {
		$$ = new conditions_def;
		$$->type = 0;
		$$->intv = 7;
		$$->left = $1;
		$$->right = $3;	
	}

%%
