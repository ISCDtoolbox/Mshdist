# Set the minimum version of cmake (and default installation rpath if needed)
cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
if (APPLE)
    set(CMAKE_MACOSX_RPATH TRUE)
endif ()
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# One may want to impose the C-compiler (for example with a specific version)
#set(CMAKE_C_COMPILER "/usr/bin/clang-3.8")
#set(CMAKE_C_COMPILER_ID Clang)

# Define some macro to set the correct date and time of compilation
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

# Define the project of the Mshdistance library depending on cmake version
if (CMAKE_VERSION VERSION_GREATER 3.0)
    project(Mshdistance LANGUAGES C)
else ()
    project(Mshdistance C)
endif ()
set(CMAKE_BUILD_TYPE "Release")

# Print a warning if the compiler is not gcc or clang (apple clang on mac)
if (NOT ${CMAKE_C_COMPILER_ID} MATCHES GNU
                 AND NOT ${CMAKE_C_COMPILER_ID} MATCHES Clang
                              AND NOT ${CMAKE_C_COMPILER_ID} MATCHES AppleClang)
    message("Warning: the Mshdistance library is currently intended to be ")
    message("compiled by gcc or clang since it has been tested merely on ")
    message("these two compilers. The commands 'cmake ..' and 'make' may not ")
    message("work.")
endif ()

# Write compilation date in a file entitled compil.date of sources/ directory
DATE(day)
TIME(hour)
file(WRITE sources/compil.date "#define COMPIL \"${day} ${hour}\"")

# Add the source files to the Mshdistance library which also uses
#    * Headers from the C standard library:
#        - assert.h
#        - ctype.h
#        - limits.h
#        - math.h (the associated static library needs to be linked)
#        - signal.h
#        - stdio.h
#        - stdlib.h
#        - string.h
#
#    * Headers from the external (dynamic) Commons library, which needs to be
#      properly linked:
#         - chrono.h
#         - libmesh5.h
#         - lplib3.h
#
#    * Headers from the local source directory:
#        - mshdist.h
#        - memory.h
#        - compil.date (created only at the compilation by cmake, cf. macro)
add_library(Mshdistance SHARED sources/boulep.c
                               sources/bucket.c
                               sources/hash.c
                               sources/inout.c
                               sources/locelt.c
                               sources/memory.c
                               sources/mshdis1_2d.c
                               sources/mshdis1_3d.c
                               #sources/mshdist.c
                               sources/queue.c
                               sources/scalem.c
                               sources/specdist_2d.c
                               sources/specdist_3d.c
                               sources/tools_2d.c
                               sources/tools_3d.c)

# Find the Commons library (all the other libraries are its dependancies)
find_library(COMMONS_LIBRARY NAMES Commons HINTS "$ENV{HOME}/lib/iscd")

# Find the Math library (libm.a on Ubuntu 16.04.4, add -lm with gcc compiler)
find_library(MATH_LIBRARY m)

# Set all the targets and properties of the Mshdistance library
target_include_directories(Mshdistance
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/sources>
        $<BUILD_INTERFACE:$ENV{HOME}/include/iscd/Commons>)
        # Remark: the user can add directly the path(s) for header file(s) in
        #         the above target_include_directories command (such as the ones
        #         associated with the Commons library)

if (CMAKE_VERSION VERSION_GREATER 3.8)
    target_compile_features(Mshdistance PRIVATE c_std_99)
elseif (CMAKE_VERSION VERSION_GREATER 3.1)
    set_property(TARGET Mshdistance PROPERTY C_STANDARD 99)
else ()
    target_compile_options(Mshdistance PRIVATE -std=gnu99)
endif ()

if (${CMAKE_C_COMPILER_ID} MATCHES GNU)
    target_compile_options(Mshdistance PRIVATE -Wno-unused-result)
    #target_compile_options(Mshdistance PRIVATE -Wpedantic -Wall -Wextra)
elseif (${CMAKE_C_COMPILER_ID} MATCHES Clang
                                  OR ${CMAKE_C_COMPILER_ID} MATCHES AppleClang)
    if (APPLE)
        # Important remark: on apple clang, the 'info' global variable of
        #                   mshdist.c causes problems when linking the
        #                   Mshdistance library
        target_link_libraries(Mshdistance PRIVATE -Wl,-U,_info)
    endif ()
#else ()
    #target_compile_options(Mshdistance PRIVATE -g)
    #target_compile_options(Mshdistance PRIVATE -Weverything)
endif ()

# Remark: the user can add directly here the path for linking the math and
#         and Commons libraries (if they have not been found with the above
#         find_library commands)
target_link_libraries(Mshdistance PUBLIC ${MATH_LIBRARY}
                                         ${COMMONS_LIBRARY})

# Set the installation properties of the Mshdistance library
# Remark: the user can change here the paths where the Mshdistance library and
#         associated header files will be installed at 'cmake install' command
install(TARGETS Mshdistance LIBRARY DESTINATION "$ENV{HOME}/lib/iscd")
install(FILES sources/memory.h
              sources/mshdist.h DESTINATION "$ENV{HOME}/include/iscd/Mshdist")

# Create and install the corresponding mshdist executable
if (CMAKE_VERSION VERSION_GREATER 3.0)
    project(mshdist LANGUAGES C)
else ()
    project(mshdist C)
endif ()
set(CMAKE_BUILD_TYPE "Release")

add_executable(mshdist sources/mshdist.c)

if (CMAKE_VERSION VERSION_GREATER 3.8)
    target_compile_features(mshdist PRIVATE c_std_99)
elseif (CMAKE_VERSION VERSION_GREATER 3.1)
    set_property(TARGET mshdist PROPERTY C_STANDARD 99)
else ()
    target_compile_options(mshdist PRIVATE -std=gnu99)
endif ()

if (${CMAKE_C_COMPILER_ID} MATCHES GNU)
    target_compile_options(mshdist PRIVATE -Wno-unused-result)
    #target_compile_options(mshdist PRIVATE -Wpedantic -Wall -Wextra)
else ()
    #target_compile_options(mshdist PRIVATE -g)
    #target_compile_options(mshdist PRIVATE -Weverything)
endif ()

target_link_libraries(mshdist PUBLIC Mshdistance)
install(TARGETS mshdist RUNTIME DESTINATION "$ENV{HOME}/bin")

