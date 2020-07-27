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
string [\'][^\']*[\']
id [A-Za-z_][A-Za-z0-9_]*
number (\-|\+)?[0-9]+
float (\-|\+)?[0-9]+\.[0-9]*
comment \-\-(\s|.)*
comment2 \/\*(\s|.)*?\*\/

%%

and    {return AND;}
as     {return AS;}
char   {return CHAR;}
create {return CREATE;}
delete {return DELETE;}
database {return DATABASE;}
drop   {return DROP;}
exit   {return EXIT;}
float  {return FLOAT;}
from   {return FROM;}
in	   {return IN;}
int	   {return INT;}
insert {return INSERT;}
into   {return INTO;}
not    {return NOT;}
null   {return NULLV;}
or     {return OR;}
set    {return SET;}
save   {return SAVE;}
select {return SELECT;}
schema {return SCHEMA;}
show   {return SHOW;}
table  {return TABLE;}
unique {return UNIQUE;}
update {return UPDATE;}
use    {return USE;}
values {return VALUES;}
where  {return WHERE;}

[.] |
[;] | 
[(] | 
[)] |
[,] |
[*] |
[=] |
[!] |
[>] |
[<] {return *yytext;}
 
{id}  {
    yylval.strval=strdup(yytext);
    for (char *pc=yylval.strval;*pc;pc++)
        *pc=tolower(*pc);
    return ID;
}
{number} {yylval.intval=atoi(yytext); return NUMBER;}
{float} {yylval.floatval=atof(yytext); return FLOATVAL;}
{ws}
{comment}
{comment2}
{string}  {
    /*remove ''*/
    yylval.strval=(char*)malloc(yyleng-1);
    memcpy(yylval.strval,yytext+1,yyleng-1);
    yylval.strval[yyleng-2]=0;
    return STRING;
}
\n     {return *yytext;}

 /*bad character*/
. {fprintf(stderr, "lex error, unexpected charactor '%s'\n", yytext);}

%%
