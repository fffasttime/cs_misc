
args = -std=c++11 -g -D_GLIBCXX_DEBUG -Wall

all : common.o simplesql.o storage.o parser.o lexer.o
	g++ common.o simplesql.o storage.o parser.o lexer.o -o simplesql $(args)
common.o : common.h common.cc
	g++ -o common.o -c common.cc $(args)
simplesql.o : simplesql.h simplesql.cc
	g++ -c simplesql.cc $(args)
storage.o : storage.h storage.cc
	g++ -c storage.cc $(args)
parser.o : simplesql.y
	bison -d -o parser.cc simplesql.y
	g++ -c parser.cc $(args)
lexer.o : simplesql.l
	flex -o lexer.cc simplesql.l
	g++ -c lexer.cc $(args)
clean : 
	rm -f lexer.* parser.* simplesql *.o
