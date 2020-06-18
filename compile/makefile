
all : parser.o lexer.o simplesql.o
	g++ -o simplesql parser.o lexer.o simplesql.o -std=c++11
simplesql.o : simplesql.h simplesql.cc
	g++ -c simplesql.cc
lexer.o : simplesql.l
	flex -o lexer.cc simplesql.l
	g++ -c lexer.cc
parser.o : simplesql.y
	bison -d -o parser.cc simplesql.y
	g++ -c parser.cc
clean : 
	rm lexer.* parser.* simplesql *.o
