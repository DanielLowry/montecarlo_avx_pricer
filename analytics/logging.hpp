#pragma once
#include <iostream>
#include <sstream>
#include <iomanip>

#ifdef ENABLE_LOGGING
#define LOG(msg) do { \
    std::stringstream ss; \
    ss << "[LOG] " << std::fixed << std::setprecision(6) << msg; \
    std::clog << ss.str() << std::endl; \
} while(0)
#else
#define LOG(msg) do {} while(0)
#endif
