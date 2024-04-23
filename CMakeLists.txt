cmake_minimum_required(VERSION 3.12)

project(Mshdist C)

# ---- Configure compiler ----------------------------------------------------
add_compile_options(
  -Wno-char-subscripts
  -Wno-unused-parameter
  -Wno-unused-variable
  -Wno-sign-compare
  -Wno-incompatible-pointer-types
  -Wno-incompatible-function-pointer-types)

# ---- ISCD::Commons ---------------------------------------------------------
add_subdirectory(
  ${PROJECT_SOURCE_DIR}/third-party/ISCD/Commons
  ${PROJECT_BINARY_DIR}/third-party/ISCD/Commons)
target_include_directories(Commons
  INTERFACE
  ${PROJECT_SOURCE_DIR}/third-party/ISCD/Commons/sources)
add_library(ISCD::Commons ALIAS Commons)

# ---- Executable ------------------------------------------------------------
file(GLOB_RECURSE Mshdist_SRCS sources/*.c)
file(GLOB_RECURSE Mshdist_HEADERS sources/*.h)

add_executable(mshdist ${Mshdist_SRCS} ${Mshdist_HEADERS})

target_link_libraries(mshdist
  PRIVATE
  m
  pthread
  ISCD::Commons)
