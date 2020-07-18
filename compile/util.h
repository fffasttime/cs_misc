#ifndef UTIL_H
#define UTIL_H

#include <string>
using std::string;

[[noreturn]] fatal(char *fmt, ...);
int printf_debug(char *fmt, ...);
int printf_info(char *fmt, ...);
int printf_error(char *fmt, ...);

enum class FieldType{
	int32,
	nchar //char(n)
};


#endif
