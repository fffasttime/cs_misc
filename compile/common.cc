#include <cstdarg>
#include <cstdio>
#include "common.h"

[[noreturn]] void fatal(const char *fmt, ...){
    va_list ap;
    va_start(ap, fmt);
    vfprintf(stderr, fmt, ap);
    va_end(ap);
    exit(1);
}

int printf_debug(const char *fmt, ...){
#ifndef NODEBUG
    va_list ap;
    va_start(ap, fmt);
    int ret=vfprintf(stdout, fmt, ap);
    va_end(ap);
    return ret;
#else
    return 0
#endif
}

int printf_info(const char *fmt, ...){
    va_list ap;
    va_start(ap, fmt);
    int ret=vfprintf(stdout, fmt, ap);
    va_end(ap);
    return ret;
}

int printf_error(const char *fmt, ...){
    va_list ap;
    va_start(ap, fmt);
    int ret=vfprintf(stderr, fmt, ap);
    va_end(ap);
    return ret;
}
