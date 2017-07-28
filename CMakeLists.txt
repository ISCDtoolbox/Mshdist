cmake_minimum_required(VERSION 2.8)
set(CMAKE_BUILD_TYPE  "Release")
INCLUDE_DIRECTORIES(  "$ENV{HOME}/include")
LINK_DIRECTORIES( "$ENV{HOME}/lib")
SET(CMAKE_CXX_FLAGS   "-w -O3")
SET(CMAKE_C_FLAGS   "-w -O3")
set(CMAKE_MACOSX_RPATH 1)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

#Advection library
project(MshDist)
file(                  GLOB_RECURSE source_files sources/*)
file(                  GLOB_RECURSE header_files sources/*.h)

FIND_LIBRARY(          M_LIB m)
find_library(          Commons NAMES Commons HINTS "$ENV{HOME}/lib")

ADD_LIBRARY(           MshDistance STATIC ${source_files})
target_link_libraries( MshDistance ${Commons} ${M_LIB})
#INSTALL(               TARGETS MshDistance LIBRARY  DESTINATION "$ENV{HOME}/lib")
INSTALL(               FILES ${header_files} DESTINATION "$ENV{HOME}/include")

#Executable
add_executable(        mshdist sources/mshdist.c)
target_link_libraries( mshdist MshDistance ${M_LIB})
INSTALL(               TARGETS mshdist RUNTIME DESTINATION "$ENV{HOME}/bin")
