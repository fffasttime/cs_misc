%{

#include "simplesql.h"

#define HINTCMD printf("SQL>")

%}

%union{
	int intval;
	char *strval;
	
}

%token AND CREATE DELETE DROP EXIT FROM SELECT TABLE WHERE
%token <intval> NUMBER
%token <strval> STRING ID INT CHAR

%left OR
%left AND

%%

statements: statements statement | statement
statement: createsql | selectsql | exitsql

selectsql: SELECT '*' FROM tables ';' '\n'{
		puts("");
		selection(NULL, $4, NULL);
		puts(""); HINTCMD;
	}
	| SELECT item_list FROM tables ';' '\n' {
		puts("");
		selection($2, $4, NULL);
		puts(""); HINTCMD;	
	}
	| SELECT '*' FROM tables WHERE conditions ';' '\n' {
		puts("");
		selection(NULL, $4, $6);
		puts("") HINTCMD;
	}
	| SELECT '*' FROM tables WHERE conditions ';' '\n' {
		puts("");
		selection($2, $4, $6);
		puts("") HINTCMD;
	}

createsql: CREATE TABLE ID '(' create_items ')' ';' '\n' {
		puts("");
		createTable($3, $5);
		puts("") HINTCMD;
	}
	/*todo: create database*/

exitsql: EXIT ';' {
		puts("");
		printf("Successfully exited\n");
		exit(0);
	}

create_item: ID INT{
		$$=(hyper_item_def *)malloc(sizeof(hyper_item_def));
		$$->field=$1;
		$$->type=0;
		$$->next=NULL;
	}
	/*todo: text field*/

create_items:  create_item { $$ = $1; }
	| create_items ',' create_item{
		$$ = $3;
		$$->next=$1;	
	}

%%
