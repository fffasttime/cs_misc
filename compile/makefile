
all : parser.cc lexer.cc 
	gcc -o simplesql parser.cc lexer.cc simplesql.cc
lexer.cc : simplesql.l
	flex -o lexer.cc simplesql.l
parser.cc : simplesql.y
	bison -o parser.cc simplesql.y
clean : 
	-rm lexer.cc parser.cc simplesql
