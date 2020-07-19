#ifndef COMMON_H
#define COMMON_H

#include <string>
#include <cassert>
using std::string;

#define DB_DATA_PATH "./data/"

[[noreturn]] void fatal(const char *fmt, ...);
int printf_debug(const char *fmt, ...);
int printf_info(const char *fmt, ...);
int printf_error(const char *fmt, ...);

enum class FieldType{
	int32,
	nchar //char(n)
};


#endif
