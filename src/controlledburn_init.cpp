// controlledburn_init.cpp
// The init function is called from R .onLoad.

#include <cpp11.hpp>

extern "C" void libgeos_init_api(void);

[[cpp11::register]]
void cpp_controlledburn_init() {
    libgeos_init_api();
}
