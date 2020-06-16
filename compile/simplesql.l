%{

#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include "y.tab.h"

%}

%option caseless

%%

and    {return AND;}
create {return CREATE;}
delete {return DELETE;}
drop   {return DROP;}
exit   {return EXIT;}
from   {return FROM;}
select {return SELECT;}
table  {return TABLE;}
where  {return WHERE;}

[;] | 
[(] | 
[)] |
[,] |
[*] |
[=] |
[!] |
[>] |
[<] |
[<=] |
[>=] |
[<>]   {return *yytext;}
 
 /*id*/
[A-Za-z][A-Za-z0-9_]*   {yylval.strval=strdup(yytext); return ID;}
 /*str*/
[\'][A-Za-z][A-Za-z0-9_]*[\']  {yylval.strval=strdup(yytext); return STRING;}
 /*num*/
[0-9]+  {yylval.intval=atoi(yytext); return NUMBER;}
\n     {return *yytext;}
[ \t]+

%%

int main(int argc, char **argv){
	printf("SQL>");	
	yyparse();
	return 0;
}

void yyerror(const char *str){
	fprintf(stderr, "Error: %s at line %d\n", str, yylineno);
}

int yywarp(){
	return 1;
}
