#pragma once
#include <iostream>

#ifdef ENABLE_LOGGING
//#error "Logging shouldn't be defined"
#define LOG(msg) do { std::clog << "[LOG] " << msg << std::endl; } while(0)
#else
#define LOG(msg) do {} while(0)
#endif
