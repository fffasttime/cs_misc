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
	string strval;
	create_item_def *crtitemval;
	select_item_def *selitemval;
	conditions_def *conval;
	table_def *tbval;
}

%error-verbose

%token AND CREATE DATABASE DELETE DROP EXIT FROM INSERT 
%token INTO SAVE SCHEMA SELECT TABLE USE WHERE
%token <intval> NUMBER
%token <strval> STRING ID INT CHAR

/* binding content */
%type <intval> comparator
%type <crtitemval> create_items create_item
%type <conval> condition conditions
%type <selitemval> item item_list
%type <tbval> tables

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
		puts("");
	}
	| SELECT item_list FROM tables {
		selection($2, $4, nullptr);
		puts("");
	}
	| SELECT '*' FROM tables WHERE conditions {
		selection(nullptr, $4, $6);
		puts("");
	}
	| SELECT item_list FROM tables WHERE conditions {
		selection($2, $4, $6);
		puts("");
	}


exitsql: EXIT {
		printf_info("info: Saving data\n");
		saveDatabase();
		printf_info("info: Successfully exited\n");
		exit(0);
	}

createsql: CREATE TABLE ID '(' create_items ')' {
		createTable($3, $5);
	}
	| CREATE DATABASE ID{
		createDatabase($3);
	}
	| CREATE SCHEMA ID{
		createDatabase($3);
	}

create_item: ID INT{
		$$ = new create_item_def;
		$$->name=$1;
		$$->type=FieldType::int32;
	}
	| CHAR '(' NUMBER ')' {
		$$->name=$1;
		$$->type=FieldType::nchar;
		$$->extra=$3;
	}

create_items:  create_item { 
		$$ = new create_item_def;
		$$.emplace_back($1);
	}
	| create_items ',' create_item{
		$$ = $1;
		$$.push_back($3);
	}

insertsql: INSERT INTO ID VALUES '(' value_list ')'{
		insertData($3, $6);
	}
	| INSERT INTO ID '(' item_list ')' VALUES '(' value_list ')'{
		insertData($3, $9, $5);
	}

value_list: value {$$ = $1;}
	| value_list ',' value{
		$$ = $3;
		$$->next = $1;
	}

value: NUMBER {
		$$ = new value_def;
		$$.emplace_back($1, FieldType::int32);
	}
	| STRING {
		$$ = new value_def;
		$$.strkey=$1;
		$$->type = FieldType::nchar;
	}

tables:	ID {
		$$ = new table_def;
		$$.emplace_back($1);
	}
	| tables ',' ID{
		$$ = $1;
		$$.emplace_back($3);			
	}


item:   ID {
		$$ = new item_def;
		$$->name = $1;
		$$->pos = nullptr;
		$$->next = nullptr;
	}

item_list: 		  item { $$ = $1; }
	| item_list ',' item { $$ = $3; }

comparator:		  '='     {$$ = 1;}
	| '>'     {$$ = 2;}
	| '<'     {$$ = 3;}
	| ">="    {$$ = 4;}
	| "<="    {$$ = 5;}
	| "<>"    {$$ = 6;}
	| '!' '=' {$$ = 6;}

condition: item comparator NUMBER {
		$$ = new conditions_def;
		$$->type = 0;
		$$->litem = $1;
		$$->intv = $3;
		$$->cmp_op = $2;
		$$->left = nullptr;
		$$->right = nullptr;
	}
	| item comparator STRING {
		$$ = new conditions_def;
		$$->type = 1;
		$$->litem = $1;
		$$->strv = $3;
		$$->cmp_op = $2;
		$$->left = nullptr;
		$$->right = nullptr;
	}

conditions: condition { $$ = $1; }
	|   '(' condition ')' { $$ = $2; }
	|   conditions AND conditions {
		$$ = new conditions_def;
		$$->cmp_op = 7;
		$$->left = $1;
		$$->right = $3;	
	}
	|   conditions OR conditions {
		$$ = new conditions_def;
		$$->cmp_op = 8;
		$$->left = $1;
		$$->right = $3;	
	}

%%
