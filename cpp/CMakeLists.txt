cmake_minimum_required(VERSION 3.16)
project(controlledburn VERSION 0.2.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

option(CONTROLLEDBURN_BUILD_TESTS "Build tests" ON)

add_library(controlledburn
  src/scanline.cpp
  src/wkb.cpp
  src/ee/box.cpp
  src/ee/coordinate.cpp
  src/ee/grid.cpp
  src/ee/measures.cpp
  src/ee/perimeter_distance.cpp
  src/ee/side.cpp
  src/ee/traversal.cpp
  src/ee/traversal_areas.cpp
)

target_include_directories(controlledburn
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src
)

# Required when the static core is linked into shared objects
# (Python extension modules, R shared libs).
set_target_properties(controlledburn PROPERTIES POSITION_INDEPENDENT_CODE ON)

target_compile_options(controlledburn PRIVATE
  $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:-Wall -Wextra>
)

if(CONTROLLEDBURN_BUILD_TESTS)
  enable_testing()
  add_executable(test_burn tests/test_burn.cpp)
  target_link_libraries(test_burn PRIVATE controlledburn)
  add_test(NAME burn COMMAND test_burn)
endif()

include(GNUInstallDirs)
install(TARGETS controlledburn EXPORT controlledburnTargets
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(DIRECTORY include/controlledburn DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
