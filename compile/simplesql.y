%{

#include "simplesql.h"

#define HINTCMD printf("SQL>")

int yylex();
int yyparse();
void yyerror(const char *str){
	printf("Error %s\n",str);
}
extern "C" int yywrap(){
	return 1;
}

int main(int argc, char **argv){
	printf("SQL>");	
	yyparse();
	return 0;
}

%}

%union{
	int intval;
	char *strval;
	create_items_def *crtitemsval;
	conditions_def *conval;
	item_def *itemval;
	table_def *tbval;
}
%start program
%token AND CREATE DELETE DROP EXIT FROM SELECT TABLE WHERE
%token <intval> NUMBER
%token <strval> STRING ID INT CHAR

/* binding content */
%type <intval> comparator
%type <crtitemsval> create_items create_item
%type <conval> condition conditions
%type <itemval> item item_list
%type <tbval> tables

%left OR
%left AND

%%

statements: statements statement | statement
statement: createsql | selectsql | exitsql

selectsql: SELECT '*' FROM tables ';' '\n'{
		puts("");
		selection(nullptr, $4, nullptr);
		puts(""); HINTCMD;
	}
	| SELECT item_list FROM tables ';' '\n' {
		puts("");
		selection($2, $4, nullptr);
		puts(""); HINTCMD;	
	}
	| SELECT '*' FROM tables WHERE conditions ';' '\n' {
		puts("");
		selection(nullptr, $4, $6);
		puts(""); HINTCMD;
	}
	| SELECT item_list FROM tables WHERE conditions ';' '\n' {
		puts("");
		selection($2, $4, $6);
		puts(""); HINTCMD;
	}

createsql: CREATE TABLE ID '(' create_items ')' ';' '\n' {
		puts("");
		createTable($3, $5);
		puts(""); HINTCMD;
	}
	/*todo: create database*/

exitsql: EXIT ';' {
		puts("");
		printf("Successfully exited\n");
		exit(0);
	}

create_item: ID INT{
		$$=(create_items_def *)malloc(sizeof(create_items_def));
		$$->field=$1;
		$$->type=0;
		$$->next=nullptr;
	}
	/*todo: text field*/

create_items:  create_item { $$ = $1; }
	| create_items ',' create_item{
		$$ = $3;
		$$->next=$1;	
	}


tables:			ID {
					$$ = new table_def;
					$$->table = $1;
					$$->next = nullptr;
				}
				| tables ',' ID{
					$$ = new table_def;
					$$->table = $3;
					$$->next = $1;				
				}


item: 			ID {
					$$ = new item_def;
					$$->field = $1;
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
