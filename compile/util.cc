#include <cstdarg>
#include <cstdio>
#incldue "util.h"

[[noreturn]] fatal(char *fmt, ...){
    va_list ap;
    va_start(ap, fmt);
    vfprintf(stderr, fmt, ap);
    va_end(ap);
    exit(1);
}

void printf_info(char *fmt, ...){
    va_list ap;
    va_start(ap, fmt);
    int ret=vfprintf(stdout, fmt, ap);
    va_end(ap);
    return ret;
}

void printf_error(char *fmt, ...){
    va_list ap;
    va_start(ap, fmt);
    int ret=vfprintf(stderr, fmt, ap);
    va_end(ap);
    return ret;
}