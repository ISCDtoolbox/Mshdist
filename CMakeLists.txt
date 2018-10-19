cmake_minimum_required(VERSION 3.5.1 FATAL_ERROR)
#set(CMAKE_C_COMPILER "/usr/bin/clang")
#set(CMAKE_C_COMPILER_ID Clang)

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
#      the Commons libraries as well as its dependencies i.e. the Math, 
#      the Threads (libpthread-stubs0-dev on Ubuntu 16.04.4, add -pthread with
#                   gcc compiler), and optionally (depending on how the Commons
#                   library was compiled)
#      the OpenMp library (libgomp1 on Ubuntu 16.04.4, add -fopenmp with gcc
#                          compiler).
#
# FIRST REMARK: when -fopenmp option is set in gcc, the preprocessor compiler
#               creates a preprocessor constant _OPENMP yyyymmdd where yyyy is
#               the year mm the month and dd the day of the preprocessor
#               compilation.
#
# SECOND REMARK: adding -lpthread option with gcc only links the pthread library
#                while the pre-defined macros do not get defined. The -pthread
#                option tells the compiler to link in the pthread library as
#                well as configure the compilation for threads. Tendancy is to
#                add -pthread option instead of -lpthread
find_library(Commons NAMES Commons HINTS "$ENV{HOME}/lib/iscd")
#find_library(M_LIB m)
#set(THREADS_PREFER_PTHREAD_FLAG ON)
#find_package(Threads REQUIRED)
#find_package(OpenMP)

# Print a warning to suggest what to do if the OpenMp library is not found
#if (NOT OPENMP_FOUND)
#    if (APPLE)
#        message("Warning: check that the libomp library has been properly ")
#        message("installed by typing for example 'brew install libomp' ")
#        message("in your terminal window.")
#        if (CMAKE_VERSION VERSION_LESS 3.12)
#            message("Warning: on mac os, your version of cmake ")
#            message("(=${CMAKE_VERSION}) cannot be currently less than ")
#            message("3.12 in order to target properly the OpenMp library.")
#            message("Please, update your version of cmake by typing for ")
#            message("example 'brew upgrade cmake' in your terminal window.")
#        endif ()
#    else ()
#        if (${CMAKE_C_COMPILER_ID} MATCHES GNU)
#            message("Warning: check that the libgomp1 library has been ")
#            message("properly installed by typing for example ")
#            message("'sudo apt-get install libgomp1' in your terminal ")
#            message("window (you need administrator rights to do that).")
#        elseif (${CMAKE_C_COMPILER_ID} MATCHES Clang)
#            message("Warning: check that the libomp-dev library has been ")
#            message("properly installed by typing for example ")
#            message("'sudo apt-get install libomp-dev' in your terminal ")
#            message("window (you need administrator rights to do that).")
#        endif ()
#    endif ()
#endif ()

# Set all the targets and properties of the Mshdistance library
target_include_directories(Mshdistance
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/sources>
        $<BUILD_INTERFACE:$ENV{HOME}/include/iscd/Commons>)
set_property(TARGET Mshdistance PROPERTY C_STANDARD 99)
target_link_libraries(Mshdistance PUBLIC ${Commons})
#                                         ${M_LIB}
#                                         Threads::Threads)
#if (OPENMP_FOUND)
#    if (CMAKE_VERSION VERSION_GREATER 3.9)
#        target_link_libraries(Mshdistance PUBLIC OpenMP::OpenMP_C)
#    else ()
#        target_link_libraries(Mshdistance PUBLIC ${OpenMP_C_FLAGS})
#        target_compile_options(Mshdistance PUBLIC ${OpenMP_C_FLAGS})
#    endif ()
#endif ()
#target_compile_options(Mshdistance PUBLIC -Wno-unused-result)
target_compile_options(Mshdistance PUBLIC -Wall -Wextra -Wpedantic -Wno-char-subscripts) 

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
#    target_compile_options(mshdist PUBLIC -Wno-unused-result) 
    target_compile_options(mshdist PUBLIC -Wall -Wextra -Wpedantic) 
endif ()
target_link_libraries(mshdist PUBLIC Mshdistance)
install(TARGETS mshdist RUNTIME DESTINATION "$ENV{HOME}/bin")

