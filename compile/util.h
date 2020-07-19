#ifndef UTIL_H
#define UTIL_H

#include <string>
#include <cassert>
using std::string;

#define DB_DATA_PATH "data/"

[[noreturn]] void fatal(char *fmt, ...);
int printf_debug(char *fmt, ...);
int printf_info(char *fmt, ...);
int printf_error(char *fmt, ...);

enum class FieldType{
	int32,
	nchar //char(n)
};


#endif
