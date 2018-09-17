/*

  Copyright (C) 2015-2018 Yuri Pirola, Simone Zaccaria

  Distributed under the MIT license.

  You should have received a copy of the MIT license along with this
  program.

*/

#ifndef _LOG_H_
#define _LOG_H_

#define LOG_LEVEL_FATAL (0)
#define LOG_LEVEL_ERROR (1)
#define LOG_LEVEL_WARN (2)
#define LOG_LEVEL_INFO (3)
#define LOG_LEVEL_DEBUG (4)
#define LOG_LEVEL_TRACE (5)
#define LOG_LEVEL_FINETRACE (6)

#ifndef LOG_THRESHOLD
#define LOG_THRESHOLD LOG_LEVEL_INFO
#endif

#ifndef LOG_PREFIX
#define LOG_PREFIX "* "
#endif


#endif // _LOG_H_


#ifdef LOG
#undef __INTERNAL_LOG
#undef LOG
#undef FATAL
#undef ERROR
#undef WARN
#undef INFO
#undef DEBUG
#undef TRACE
#undef FINETRACE
#endif

#ifdef LOG_MSG

#include <string>
#include <iomanip>
#include <iostream>

#define MAX_LEN_FUNC_NAME 12
#define MAX_LEN_FILE_NAME 12

#define __LOG_PREFIXES__LOG_LEVEL_FATAL "FATAL"
#define __LOG_PREFIXES__LOG_LEVEL_ERROR "ERROR"
#define __LOG_PREFIXES__LOG_LEVEL_WARN  "WARN "
#define __LOG_PREFIXES__LOG_LEVEL_INFO  "INFO "
#define __LOG_PREFIXES__LOG_LEVEL_DEBUG "DEBUG"
#define __LOG_PREFIXES__LOG_LEVEL_TRACE "TRACE"
#define __LOG_PREFIXES__LOG_LEVEL_FINETRACE "FTRAC"

#define LOG(level, ...) __INTERNAL_LOG(level, __LOG_PREFIXES__ ## level, __VA_ARGS__, "")

#define ALWAYS_LOG(level, ...) __INTERNAL_ALWAYS_LOG(__LOG_PREFIXES__ ## level, __VA_ARGS__, "")

#define __INTERNAL_LOG(level, prefix, format, ...) do {                 \
        if (level<=LOG_THRESHOLD) {                                     \
            __INTERNAL_ALWAYS_LOG(prefix, format, __VA_ARGS__);         \
        }                                                               \
    } while (0)

#define __INTERNAL_ALWAYS_LOG(prefix, format, ...) do {                 \
        time_t ttNow = time(NULL); \
        tm tmNow = *localtime(&ttNow);                                  \
        std::string __str_func__=std::string(__func__);                 \
        std::string __my_internal_funz__=__str_func__.substr(0,MAX_LEN_FUNC_NAME); \
        std::string __str_FILE__=std::string(__FILE__);                 \
        std::cerr << LOG_PREFIX << prefix << "(";                       \
        std::cerr << std::setw(MAX_LEN_FUNC_NAME) << std::setfill(' '); \
        std::cerr << std::left << __my_internal_funz__;                 \
        std::cerr << ":";                                               \
        std::cerr << std::setw(MAX_LEN_FILE_NAME) << std::setfill(' '); \
        std::cerr << ((__str_FILE__.length() > MAX_LEN_FILE_NAME) ?     \
          __str_FILE__.substr(__str_FILE__.length()-MAX_LEN_FILE_NAME) : \
          __str_FILE__);                                                 \
        std::cerr << ":" << std::left << std::setw(4) << __LINE__ << std::right << ") "; \
        std::cerr << std::setfill('0') << std::setw(2) << tmNow.tm_hour << ":"; \
        std::cerr << std::setfill('0') << std::setw(2) << tmNow.tm_min << ":";  \
        std::cerr << std::setfill('0') << std::setw(2) << tmNow.tm_sec << " | ";\
        std::cerr << format  << __VA_ARGS__ << std::endl;               \
    } while(0)
#else

#define LOG(level, prefix, ...) do { } while (0)
#define ALWAYS_LOG(level, prefix, ...) do { } while (0)

#endif


#if (LOG_LEVEL_FATAL <= LOG_THRESHOLD) && defined LOG_MSG
#define LOG_FATAL_ENABLED
#define FATAL(...) LOG(LOG_LEVEL_FATAL, __VA_ARGS__)
#else
#undef LOG_FATAL_ENABLED
#define FATAL(...) do { } while (0)
#endif

#if (LOG_LEVEL_ERROR <= LOG_THRESHOLD) && defined LOG_MSG
#define LOG_ERROR_ENABLED
#define ERROR(...) LOG(LOG_LEVEL_ERROR, __VA_ARGS__)
#else
#undef LOG_ERROR_ENABLED
#define ERROR(...) do { } while (0)
#endif

#if (LOG_LEVEL_WARN <= LOG_THRESHOLD) && defined LOG_MSG
#define LOG_WARN_ENABLED
#define WARN(...) LOG(LOG_LEVEL_WARN, __VA_ARGS__)
#else
#undef LOG_WARN_ENABLED
#define WARN(...) do { } while (0)
#endif

#if (LOG_LEVEL_INFO <= LOG_THRESHOLD) && defined LOG_MSG
#define LOG_INFO_ENABLED
#define INFO(...) LOG(LOG_LEVEL_INFO, __VA_ARGS__)
#else
#undef LOG_INFO_ENABLED
#define INFO(...) do { } while (0)
#endif

#if (LOG_LEVEL_DEBUG <= LOG_THRESHOLD) && defined LOG_MSG
#define LOG_DEBUG_ENABLED
#define DEBUG(...) LOG(LOG_LEVEL_DEBUG, __VA_ARGS__)
#else
#undef LOG_DEBUG_ENABLED
#define DEBUG(...) do { } while (0)
#endif

#if (LOG_LEVEL_TRACE <= LOG_THRESHOLD) && defined LOG_MSG
#define LOG_TRACE_ENABLED
#define TRACE(...) LOG(LOG_LEVEL_TRACE, __VA_ARGS__)
#else
#undef LOG_TRACE_ENABLED
#define TRACE(...) do { } while (0)
#endif

#if (LOG_LEVEL_FINETRACE <= LOG_THRESHOLD) && defined LOG_MSG
#define LOG_FINETRACE_ENABLED
#define FINETRACE(...) LOG(LOG_LEVEL_FINETRACE, __VA_ARGS__)
#else
#undef LOG_FINETRACE_ENABLED
#define FINETRACE(...) do { } while (0)
#endif
