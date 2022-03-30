cmake_minimum_required(VERSION 3.12)
project(Mshdist)

file(GLOB_RECURSE sources sources/*.c)
file(GLOB_RECURSE headers sources/*.h)

add_library(libMshdist ${sources})

find_library(Commons
  NAMES
  Commons
  HINTS
  ${ISCD_COMMONS_LIBRARY}
  )

find_library(M_LIB m)

target_include_directories(libMshdist PUBLIC ${ISCD_COMMONS_INCLUDE_DIR})
target_link_libraries(libMshdist
  PUBLIC
  ${Commons}
  pthread
  ${M_LIB}
  )

#Executable
add_executable(        mshdist sources/mshdist.c)
target_link_libraries( mshdist libMshdist ${Commons})

# This install is really borked, so I do not want to use it
# install(               TARGETS mshdist RUNTIME  DESTINATION "bin")
# install(               TARGETS libMshdist LIBRARY DESTINATION "lib")
# install(               FILES   ${headers} DESTINATION "include")
# set_target_properties( mshdist PROPERTIES INSTALL_RPATH "/usr/local/lib")

