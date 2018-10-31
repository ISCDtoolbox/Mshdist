# mshdist [![](https://travis-ci.org/ISCDtoolbox/Mshdist.svg?branch=test_future_update)](https://travis-ci.org/ISCDtoolbox/Mshdist)

The mshdist software is a simple algorithm written in C that is able to generate the signed distance function of given objects in two and three space dimensions.

## Install (Linux & Mac OS)

First, you will need to install the ISCD [Commons Library](https://github.com/ISCDtoolbox/Commons) on your system. Please refer to the instructions provided on the ISCD Commons Library page in order to properly install this library.

Then, in a command prompt, first navigate to the directory where you want to save the files (thanks to the 'cd' and 'ls' commands) and type:
```
git clone https://github.com/ISCDtoolbox/Mshdist.git
```
to clone the repository on your computer (we recall that it is also possible to download the files manually from github in a *.zip archive format if the git software is not installed on your computer, just go to the upper right part of the current page, click on the button Clone or download, and then Download ZIP). Finally, type successively in the command prompt
```
cd Mshdist
```
to navigate to the downloaded Mshdist/ directory;
```
mkdir build
```
to create the build/ directory
```
cd build
```
to enter the build/ directory
```
cmake ..
```
to generate the Makefile with cmake
```
make
```
to compile the project using make
```
make install
```
to install the mshdist software on your computer. Depending on where are located your installation paths, you may also need the administrator rights in order to do that. In this case, you need to type in the command prompt
```
sudo make install
```

### Installation paths
You may have installed the Commons library at a different location than the default one (lib/iscd subdirectory of your home directory). In this case you need to modify the CMakeLists.txt file located in the Mshdist/ directory, and correct the line
```
find_library(COMMONS_LIBRARY NAMES Commons HINTS "$ENV{HOME}/lib/iscd")
```
by providing the right path for searching the Commons library. Similarly, you may have also modified the installation path for the associated public header files (by default in the include/iscd/Commons/ subdirectory of your home directory). In this case, you also have to change the default location in the follwing lines of the CMakeLists.txt file
```
target_include_directories(Mshdistance
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/sources>
        $<BUILD_INTERFACE:$ENV{HOME}/include/iscd/Commons>)
```
The mshdist software is compiled by first creating a static Mshdistance library which is then linked to it. The libMshdistance.a is installed in the lib/iscd subdirectory of your home directory by default but you can again specify another location by modifying the following line of the CMakeLists.txt file
```
install(TARGETS Mshdistance ARCHIVE DESTINATION "$ENV{HOME}/lib/iscd")
```
Similarly, the associated public header files are saved by default in the include/iscd/Mshdist/ subdirectory of your home directory, which can be changed by modifying the lines
```
install(FILES sources/memory.h
              sources/mshdist.h DESTINATION "$ENV{HOME}/include/iscd/Mshdist")
```
Finally, the mshdist software is installed by default in the bin/ subdirectory of your home directory but this can be changed by modyfing the path at the line
```
install(TARGETS mshdist RUNTIME DESTINATION "$ENV{HOME}/bin")
``` 

#### Usage

* A short documentation is supplied (see the [documentation](documentation) folder).

#### Quickstart

*

#### Authors & contributors
* mshdist is developped and maintained by Charles Dapogny (Université Joseph Fourier) and Pascal Frey (Université Pierre et Marie Curie).
* Contributors to this project are warmly welcomed.

#### License
* mshdist is given under the [terms of the GNU Lesser General Public License] (LICENSE.md).

* If you use mshdist in your work, please refer to the journal article:

C. Dapogny, P. Frey, _Computation of the signed distance function to a discrete contour on adapted triangulation_, Calcolo, Volume 49, Issue 3, pp.193-219 (2012).
