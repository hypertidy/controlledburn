// controlledburn_init.cpp
// Placeholder init function called from R .onLoad.

#include <cpp11.hpp>

[[cpp11::register]]
void cpp_controlledburn_init() {
    // Previously initialized libgeos; no longer needed since
    // burn_sparse (the only GEOS user) has been removed.
}
