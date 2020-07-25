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
	initDB();
	hintCMD();
	yyparse();
	return 0;
}

%}

%union{
	int intval;
	char *strval;
	create_item_def *crtitemval;
	create_item_def_unit crtitemvalu;
	select_item_def *selitemval;
	value_def *valdef;
	value_def_unit valdefu;
	conditions_def *conval;
	table_def *tbval;
}

%error-verbose

%token AND CREATE CHAR DATABASE DELETE DROP EXIT FROM INSERT 
%token INT INTO SAVE SCHEMA SELECT TABLE USE VALUES WHERE
%token <intval> NUMBER
%token <strval> STRING ID

/* binding content */
%type <intval> comparator
%type <crtitemvalu> create_item
%type <crtitemval> create_items
%type <conval> condition conditions condition_item
%type <selitemval> items
%type <tbval> tables
%type <valdefu> value
%type <valdef> value_list

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

statement: createsql | selectsql | exitsql | insertsql | %empty 
	| usestmt | savestmt

usestmt: USE ID {
		useDatabase($2);
	}
savestmt: SAVE{
		saveDatabase();
	}

selectsql: SELECT '*' FROM tables {
		selection(nullptr, $4, nullptr);
	}
	| SELECT items FROM tables {
		selection($2, $4, nullptr);
	}
	| SELECT '*' FROM tables WHERE conditions {
		selection(nullptr, $4, $6);
	}
	| SELECT items FROM tables WHERE conditions {
		selection($2, $4, $6);
	}


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

create_item: ID INT {
		$$.name=$1;
		$$.type=FieldType::int32;
	}
	| ID CHAR '(' NUMBER ')' {
		$$.name=$1;
		$$.type=FieldType::nchar;
		$$.extra=$4;
	}

create_items:  create_item { 
		$$ = new create_item_def({$1});
	}
	| create_items ',' create_item{
		$$ = $1;
		$$->emplace_back($3);
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

value: NUMBER {
		$$.value.intval = $1;
		$$.type = FieldType::int32;
	}
	| STRING {
		$$.value.strval = $1;
		$$.type = FieldType::nchar;
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

items: ID { $$ = new select_item_def({$1}); }
	| items ',' ID { 
		$$ = $1;
		$$->emplace_back($3);
	}

comparator:		  '='     {$$ = 1;}
	| '>'     {$$ = 2;}
	| '<'     {$$ = 3;}
	| ">="    {$$ = 4;}
	| "<="    {$$ = 5;}
	| "<>"    {$$ = 6;}
	| '!' '=' {$$ = 6;}

condition_item: ID {
		$$ = new conditions_def;
		$$->type = 3;
		$$->strv = $1;
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

condition: condition_item { $$ = $1; }
	| condition comparator condition_item {
		$$ = new conditions_def;
		$$->type = 0;
		$$->intv = $2;
		$$->left = $1;
		$$->right = $3;	
	}

conditions: condition { $$ = $1; }
	|   '(' conditions ')' { $$ = $2; }
	|   conditions AND condition {
		$$ = new conditions_def;
		$$->type = 0;
		$$->intv = 7;
		$$->left = $1;
		$$->right = $3;	
	}
	|   conditions OR condition {
		$$ = new conditions_def;
		$$->type = 0;
		$$->intv = 8;
		$$->left = $1;
		$$->right = $3;	
	}

%%
