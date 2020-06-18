
all : parser.o lexer.o simplesql.o
	g++ -o simplesql parser.o lexer.o simplesql.o -std=c++11
simplesql.o : simplesql.h simplesql.cc
	g++ -c simplesql.cc -std=c++11
lexer.o : simplesql.l
	flex -o lexer.cc simplesql.l
	g++ -c lexer.cc -std=c++11
parser.o : simplesql.y
	bison -d -o parser.cc simplesql.y
	g++ -c parser.cc -std=c++11
clean : 
	rm lexer.* parser.* simplesql *.o
