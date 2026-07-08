cmake_minimum_required(VERSION 3.16)
project(controlledburn_py LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(pybind11 CONFIG REQUIRED)

# The C++ core lives at the repository root under cpp/. This works for
# git installs and local development; a PyPI sdist would need the core
# vendored into the python/ tree (a sync step like tools/sync-core.sh
# on the R side) -- deliberate non-goal until first release.
set(CONTROLLEDBURN_BUILD_TESTS OFF CACHE BOOL "" FORCE)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../cpp controlledburn_core)

pybind11_add_module(_core src/_core.cpp)
target_link_libraries(_core PRIVATE controlledburn)

install(TARGETS _core DESTINATION controlledburn)
