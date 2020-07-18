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
	db_init();
	hintCMD();
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

%error-verbose

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

statements: statements statementline | statementline
statementline: statement ';' 
			| statement '\n' {hintCMD();} 
			| error '\n' { /*error recovery*/
				yyerrok;
				hintCMD();
			}

statement: createsql | selectsql | exitsql | %empty

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
		printf("Successfully exited\n");
		exit(0);
	}

createsql: CREATE TABLE ID '(' create_items ')' {
		createTable($3, $5);
		puts("");
	}
	/*todo: create database*/
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
