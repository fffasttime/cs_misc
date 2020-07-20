%{

#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <ctype.h>
#include "simplesql.h"
#include "parser.hh"

%}

%option caseless
/*%option c++*/

ws [ \t]+
string [\'][A-Za-z][A-Za-z0-9_]*[\']
id [A-Za-z][A-Za-z0-9_]*
number [0-9]+
comment \-\-[^\n]*
comment2 \/\*.*\*\/

%%

and    {return AND;}
char   {return CHAR;}
create {return CREATE;}
delete {return DELETE;}
database {return DATABASE;}
drop   {return DROP;}
exit   {return EXIT;}
from   {return FROM;}
int	   {return INT;}
insert {return INSERT;}
into   {return INTO;}
save   {return SAVE;}
select {return SELECT;}
schema {return SCHEMA;}
table  {return TABLE;}
use    {return USE;}
values {return VALUES;}
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
"<=" |
">=" |
"<>"   {return *yytext;}
 
{id}  {
    yylval.strval=strdup(yytext); ; 
    for (char *pc=yylval.strval;*pc;pc++)
        *pc=tolower(*pc);
    return ID;
}
{string}  {yylval.strval=strdup(yytext); return STRING;}
{number} {yylval.intval=atoi(yytext); return NUMBER;}
{ws}
{comment}
{comment2}
\n     {return *yytext;}

 /*bad character*/
. {fprintf(stderr, "lex error, unexpected charactor '%s'\n", yytext);}

%%
