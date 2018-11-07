# mshdist [![](https://travis-ci.org/ISCDtoolbox/Mshdist.svg?branch=test_future_update)](https://travis-ci.org/ISCDtoolbox/Mshdist)

The mshdist software is a simple algorithm written in C that is able to
generate the signed distance function of given objects in two and three space
dimensions.

## Install (Linux & Mac OS)

First, you will need to install the ISCD
[Commons Library](https://github.com/ISCDtoolbox/Commons) on your system.
Please refer to the instructions provided on the ISCD Commons Library page in
order to properly install this required library.

Then, in a command prompt, first navigate to the directory where you want to
save the files (thanks to the `cd` and `ls` commands) and type:
```
git clone https://github.com/ISCDtoolbox/Mshdist.git
```
to clone the repository on your computer (we recall that it is also possible to
download the files manually from github in a *\*.zip* archive format if `git`
is not installed, just go to the upper right part of the current
[project page](https://github.com/ISCDtoolbox/Mshdist), click on the button
`Clone or download`, and then *Download ZIP*). Finally, type successively in
the command prompt
```
cd Mshdist
```
to navigate to the downloaded *Mshdist/* directory;
```
mkdir build
```
to create the *build/* directory
```
cd build
```
to enter the *build/* directory
```
cmake ..
```
to generate the *Makefile* with `cmake`
```
make
```
to compile the project using `make`
```
make install
```
to install the mshdist software on your computer. Depending on where are
located your installation paths, you may also need the administrator rights in
order to do that. In this case, you need to type in the command prompt
```
sudo make install
```
If an error occurs during the installation process, it is recommended to delete
the content of the *build/* directory and start again the `cmake ..` and `make`
commands, which only generated files in this local *build/* directory.

## Installation paths
You may have installed the Commons library at a different location than the
default one (*lib/iscd* subdirectory of your *HOME* directory). In this case
you need to modify the *CMakeLists.txt* file located in the *Mshdist/*
directory, and correct the line
```
find_library(COMMONS_LIBRARY NAMES Commons HINTS "$ENV{HOME}/lib/iscd")
```
by providing the right path for searching the Commons library. Similarly, you
may have also modified the installation path for the associated public header
files (by default in the *include/iscd/Commons/* subdirectory of your *HOME*
directory). In this case, you also have to change the default location in the
following lines of the *CMakeLists.txt* file
```
target_include_directories(Mshdistance
    PUBLIC
        $<BUILD_INTERFACE:$ENV{HOME}/include/iscd/Commons>
```
The mshdist software is compiled by first creating a static Mshdistance library
which is then linked to it. The *libMshdistance.a* is installed in the
*lib/iscd* subdirectory of your *HOME* directory by default but you can again
specify another location by modifying the following line of the *CMakeLists.txt*
file
```
install(TARGETS Mshdistance ARCHIVE DESTINATION "$ENV{HOME}/lib/iscd")
```
Similarly, the associated public header files are saved by default in the
*include/iscd/Mshdist/* subdirectory of your *HOME* directory, which can be
changed by modifying the lines
```
install(FILES sources/memory.h
              sources/mshdist.h DESTINATION "$ENV{HOME}/include/iscd/Mshdist")
```
Finally, the mshdist software is installed by default in the *bin/* subdirectory
of your *HOME* directory but this can be changed by modifying the path at the line
```
install(TARGETS mshdist RUNTIME DESTINATION "$ENV{HOME}/bin")
```

## Usage

Once the mshdist software installed, you can print the usage by typing in a
command prompt
```
mshdist -?
```
A short
[documentation](https://github.com/ISCDtoolbox/Mshdist/blob/master/documentation/mshdistguide.pdf)
is also supplied with examples (see the
[documentation](https://github.com/ISCDtoolbox/Mshdist/tree/master/documentation)
folder). We provide below the typical output of the above usage command.
```
mshdist -?
  -- MSHDIST, Release 1.1b (June 21, 2010)
     Copyright (c) LJLL, 2010
     2018-10-31 16:26:19
usage: mshdist [-v[n]] [-h] file1[.mesh] [file2[.mesh]] options

** Generic options :
-d      Turn on debug mode
-?      Print this message
-dt     Time stepping (hmin)
-it n   Max number of iterations
-ncpu n Use n CPUs
-r res  Residual
-v [n]  Tune level of verbosity

  ELAPSED TIME  0.000s
```

## Examples

We have provided an example for each different mode that can be treated by the
mshdist software. In a command prompt, first navigate to the
*documentation/Examples/* subdirectory of the *Mshdist/* directory in order to
test the examples as follows.

The *\*.mesh* and their associated *\*.sol* files can be vizualized by
installing the [medit software](https://github.com/ISCDtoolbox/Medit). In this
case, just type in a command prompt
```
medit nameOfTheMeshIWantToPlot.mesh
```
There is a French
[inline HTML documentation](https://www.ljll.math.upmc.fr/frey/logiciels/Docmedit.dir/index.html)
written by [C. Dobrzynski](https://www.math.u-bordeaux.fr/~dobrzyns/)
(Université de Bordeaux) that describes how to use medit, and also a
[technical report](https://www.ljll.math.upmc.fr/frey/publications/RT-0253.pdf)
written by [Pascal Frey](https://www.ljll.math.upmc.fr/frey/)
(Sorbonne Université).

For example, you can type the `F1` or `F2` key to start or end a plane cut, and
the `m` shortcut plots (if available) the *\*.sol* values associated with the
*\*.mesh* file. You can type `h` during the execution of medit in order to have
a summary of all it main features. You can also type in a command prompt
```
medit -h
```
to print the usage of medit. More generally, we refer to the
[github page](https://github.com/ISCDtoolbox/Medit) for further details about
the medit software.

### In two space dimensions

#### With two *.mesh files

You can generate the signed distance function of the French border (edges in
*frmap.mesh* file) at the vertices of a square (given in *carre.mesh* file) by
typing in the command prompt
```
mshdist carre.mesh frmap.mesh
```
On success, the values of the signed distance function at the vertices of the
mesh will be saved in an output *carre.sol* file. Here is an example of output
of the mshdist software
```
mshdist carre.mesh frmap.mesh
  -- MSHDIST, Release 1.1b (June 21, 2010)
     Copyright (c) LJLL, 2010
     2018-10-31 16:26:19

  -- INPUT DATA
  %% carre.mesh OPENED
  -- READING DATA FILE carre.mesh
  %% frmap.mesh OPENED
  -- READING DATA FILE frmap.mesh
     NUMBER OF GIVEN VERTICES     107990     18138
     NUMBER OF GIVEN ELEMENTS     214778      1796
  -- DATA READING COMPLETED.     0.436s

  &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
   MODULE MSHDIST-LJLL : 1.1b (June 21, 2010)
  &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

  -- PHASE 1 : ANALYSIS
  -- PHASE 1 COMPLETED.     0.023s

  -- PHASE 2 : DISTANCING
  ** Initialization
     distance
     323 correction(s)
  ** Sign identification
     4 connected component(s)
  %% 66 elements with a vertex uninitialized
  %% 12 corrected vertices
  ** Propagation [1 cpu]
     Residual 2.782963E-07 after 174 iterations
  -- PHASE 2 COMPLETED.     4.342s

  &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
   END OF MODULE MSHDIST
  &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

  -- WRITING DATA FILE carre
     dd : 1.000000 (scaling factor used to unscale solutions)
  %% carre.sol OPENED
  -- WRITING COMPLETED

  ELAPSED TIME  4.864s
```

#### With a single *.mesh file and an initial associated *.sol file

Starting from a *bat.sol* file (saved as *batInit.sol* since *bat.sol* will be
overwritten) that gives a level-set function of a *bat.mesh* file (negative
inside, positive outside the domain, zero on the boundary), you can modify the
level-set function so that it become the signed distance function associated
with the bat domain (i.e. a level-set function with a unitary gradient norm).

Since such an *\*.sol* is overwritten by mshdist during the renormalization
process, we have saved the initial *\*.sol* file as *batInit.sol*. Hence, in
order to perform this test, type successively in a command prompt
```
cp batInit.sol bat.sol
mshdist bat.mesh
```
On success, the *bat.sol* file will be overwritten as the signed distance
function of the initial bat domain (that was implicitly encoded in the
*batInit.sol* file)

#### With a single *.mesh already encoding an internal domain

You can generate the signed distance function associated with a domain that is
already contained in a *thks.mesh* file (exterior/internal triangles are
labelled 2/3, edges of the internal boundary of the word thanks are labelled 5)
by typing in the command prompt
```
mshdist thks.mesh -dom
```

### In three space dimensions

You recall that you can decrease the computational time by using the
multi-threading offered by the `-ncpu` flag option of mshdist. You can obtain
the total number of cpu of your system by typing in a command prompt
```
lscpu
```
on Linux or
```
sysctl hw
```
on Mac OS systems.

#### With two *.mesh files

You can generate the signed distance function of an egg-like shape (triangles
in *mpd3d.mesh* file) at the vertices of a cube (given in *cube3d.mesh* file)
by typing in the command prompt
```
mshdist cube3d.mesh mpd3d.mesh -noscale -ncpu 1
```
On success, the values of the signed distance function at the vertices of the
mesh will be saved in an output *cube3d.sol file*. We refer to the
[documentation](https://github.com/ISCDtoolbox/Mshdist/blob/master/documentation/mshdistguide.pdf)
for a description of the `-noscale` option flag.

#### With a single *.mesh file and an initial associated *.sol file

Starting from a *ball3d.sol* file (saved as *ball3dInit.sol* since *ball3d.sol*
will be overwritten) that gives a level-set function of a *ball3d.mesh* file
(negative inside, positive outside the domain, zero on the boundary), you can
modify the level-set function so that it become the signed distance function
associated with the egg-shaped domain (i.e. a level-set function with a unitary
gradient norm).

Since such an *\*.sol* is overwritten by mshdist during the renormalization
process, we have saved the initial *\*.sol* file as *ball3dInit.sol*. Hence, in
order to perform this test, type successively in a command prompt
```
cp ball3dInit.sol ball3d.sol
mshdist ball3d.mesh -ncpu 1
```
On success, the *ball3d.sol* file will be overwritten as the signed distance
function of the initial egg-shaped domain (that was implicitly encoded in the
*ball3dInit.sol* file)

#### With a single *.mesh already encoding an internal domain

You can generate the signed distance function associated with a domain that is
already contained in a *saddle3d.mesh* file (exterior/internal triangles are
labelled 2/3, edges of the internal boundary of the saddle are labelled 10) by
typing in the command prompt
```
mshdist saddle3d.mesh -dom -ncpu 1
```

## Authors & contributors

The mshdist software is developped and maintained by
[Charles Dapogny](http://www-ljk.imag.fr/membres/Charles.Dapogny/)
(Université Joseph Fourier) and
[Pascal Frey](https://www.ljll.math.upmc.fr/frey/) (Sorbonne Université).

The associated github page has also been maintained by
[Loïc Norgeot](https://fr.linkedin.com/in/lnorgeot) and
[Jérémy Dalphin](http://pelikano.free.fr/JeremyDalphin.html).

Contributors to this project are warmly welcomed.

## Credits

The github page of the

* [Commons library](https://github.com/ISCDtoolbox/Commons)
* [Medit software](https://github.com/ISCDtoolbox/Medit)

and the link to the

* [Journal article](https://link.springer.com/article/10.1007/s10092-011-0051-z)

describing the mshdist underlying algorithm.

## License

* The mshdist software is given under the
[GNU Lesser General Public License](https://github.com/ISCDtoolbox/Mshdist/blob/master/LICENSE).

* If you use mshdist in your work, please refer to the
[journal article](https://link.springer.com/article/10.1007/s10092-011-0051-z)
as follows:

C. Dapogny, P. Frey,
_Computation of the signed distance function to a discrete contour on adapted triangulation_,
Calcolo, Volume 49, Issue 3, pp.193-219 (2012).

