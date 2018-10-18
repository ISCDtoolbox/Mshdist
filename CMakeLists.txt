cmake_minimum_required(VERSION 3.5.1 FATAL_ERROR)

# Set the installation rpath
if (APPLE)
    set(CMAKE_MACOSX_RPATH TRUE)
endif ()
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# Define some macro to set the correct date and time at compilation
macro (DATE RESULT)
  execute_process(COMMAND "date" "+%Y-%m-%d" OUTPUT_VARIABLE ${RESULT}
                                               OUTPUT_STRIP_TRAILING_WHITESPACE)
endmacro ()

macro (TIME RESULT)
  execute_process(COMMAND "date" "+%H:%M:%S" OUTPUT_VARIABLE ${RESULT}
                                               OUTPUT_STRIP_TRAILING_WHITESPACE)
endmacro ()

# Print a warning if the os is not unix (linux and mac to be more precise)
if (NOT UNIX)
    message("Warning: the Mshdistance library is currently intended for Unix ")
    message("systems since it has been tested merely on Linux and Mac OS. ")
    message("The commands 'cmake ..' and 'make' may not work.")
endif ()

# Print a warning if the compiler is not gcc or clang (apple clang on mac)
if (NOT ${CMAKE_C_COMPILER_ID} MATCHES GNU
                 AND NOT ${CMAKE_C_COMPILER_ID} MATCHES Clang
                              AND NOT ${CMAKE_C_COMPILER_ID} MATCHES AppleClang)
    message("Warning: the Mshdistance library is currently intended to be ")
    message("compiled by gcc or clang since it has been tested merely on ")
    message("these two compilers. The commands 'cmake ..' and 'make' may not ")
    message("work.")
endif ()

project(Mshdistance LANGUAGES C)
set(CMAKE_BUILD_TYPE "Release")

# Write compilation date in a file entitled compil.date of sources/ directory
DATE(day)
TIME(hour)
file(WRITE sources/compil.date "#define COMPIL \"${day} ${hour}\"")

# Add the source files to the Mshdistance library
add_library(Mshdistance SHARED sources/boulep.c
                               sources/bucket.c
                               sources/hash.c
                               sources/inout.c
                               sources/locelt.c
                               sources/memory.c
                               sources/mshdis1_2d.c
                               sources/mshdis1_3d.c
                               sources/queue.c
                               sources/scalem.c
                               sources/specdist_2d.c
                               sources/specdist_3d.c
                               sources/tools_2d.c
                               sources/tools_3d.c)

# The source files of the Mshdistance library uses:
#
#    - Headers from the C standard library: assert.h
#                                           ctype.h
#                                           limits.h
#                                           math.h (it needs to be statically
#                                                                        linked)
#                                           signal.h
#                                           stdio.h
#                                           stdlib.h
#                                           string.h
#
#    - Headers from the external Commons library, which needs to be dynamically
#      linked: chrono.h
#              libmesh5.h
#              lplib3.h
#
#    - Headers from the local sources directory: mshdist.h
#                                                memory.h
#                                                compil.date (created only at
#                                                             cmake compilation)
#
# Find the Math (libm.a on Ubuntu 16.04.4, add -lm with gcc compiler) and the
#          Commons libraries
find_library(M_LIB m)
find_library(Commons NAMES Commons HINTS "$ENV{HOME}/lib/iscd")

# Set all the targets and properties of the Mshdistance library
target_include_directories(Mshdistance
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/sources>
        $<BUILD_INTERFACE:$ENV{HOME}/include/iscd/Commons>)
set_property(TARGET Mshdistance PROPERTY C_STANDARD 99)
target_link_libraries(Mshdistance PUBLIC ${M_LIB}
                                         ${Commons})

# Set the installation properties of the Mshdistance library
install(TARGETS Mshdistance LIBRARY DESTINATION "$ENV{HOME}/lib/iscd")
install(FILES sources/memory.h
              sources/mshdist.h DESTINATION "$ENV{HOME}/include/iscd/Mshdist")

# Create and install the corresponding mshdist executable
project(mshdist LANGUAGES C)
set(CMAKE_BUILD_TYPE "Release")
add_executable(mshdist sources/mshdist.c)
set_property(TARGET mshdist PROPERTY C_STANDARD 99)
if (${CMAKE_C_COMPILER_ID} MATCHES GNU)
    target_compile_options(mshdist PUBLIC -Wno-unused-result) 
endif ()
target_link_libraries(mshdist PUBLIC Mshdistance)
install(TARGETS mshdist RUNTIME DESTINATION "$ENV{HOME}/bin")

