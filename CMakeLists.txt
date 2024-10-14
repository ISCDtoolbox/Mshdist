cmake_minimum_required(VERSION 2.8)
set(CMAKE_BUILD_TYPE  "Release")
SET(CMAKE_C_FLAGS   "-w -O3")
set(CMAKE_MACOSX_RPATH 1)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

project(libMshdist)
file(                  GLOB_RECURSE sources sources/*)
file(                  GLOB_RECURSE headers sources/*.h)
add_library(           libMshdist SHARED ${sources})
set_target_properties( libMshdist PROPERTIES PREFIX "")

find_library(          Commons NAMES Commons HINTS "$ENV{HOME}/lib" required)
FIND_LIBRARY(          M_LIB m)
SET(                   EXTRA_LIBS ${M_LIB} "-lpthread")
include_directories(   "$ENV{HOME}/include")
target_link_libraries( libMshdist ${EXTRA_LIBS} ${Commons})

#Executable
add_executable(        mshdist sources/mshdist.c)
target_link_libraries( mshdist libMshdist ${Commons})

install(               TARGETS mshdist RUNTIME  DESTINATION "bin")
install(               TARGETS libMshdist LIBRARY DESTINATION "lib")
install(               FILES   ${headers} DESTINATION "include")
set_target_properties( mshdist PROPERTIES INSTALL_RPATH "/usr/local/lib")

enable_testing()
include(CTest)

set( MSHDIST_TEST_SOURCE_DIR "${CMAKE_SOURCE_DIR}/documentation/Examples" )

ADD_TEST(NAME frmap
  COMMAND $<TARGET_FILE:mshdist>
  ${MSHDIST_TEST_SOURCE_DIR}/france/carre.mesh
  ${MSHDIST_TEST_SOURCE_DIR}/france/frmap.mesh
  )

ADD_TEST(NAME handle
  COMMAND $<TARGET_FILE:mshdist>
  ${MSHDIST_TEST_SOURCE_DIR}/handle/cube.mesh
  ${MSHDIST_TEST_SOURCE_DIR}/handle/2handle.mesh
  -ncpu 2 -it 100
  )
