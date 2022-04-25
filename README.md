# mshdist [![continuous integration](https://github.com/ISCDtoolbox/Mshdist/actions/workflows/ci.yml/badge.svg)](https://github.com/ISCDtoolbox/Mshdist/actions/workflows/ci.yml)
mshdist is a simple algorithm to generate the signed distance function to given objects in two and three space dimensions.

#### Installation

1. you will need to install the [ISCD Commons Library](https://github.com/ISCDtoolbox/Commons) on your system.
Please refer to the instructions provided on the ISCD Commons Library page in order to install this library.

2. download the zip archive of NavierStokes or clone this repository:

   ` git clone https://github.com/ISCDtoolbox/Mshdist.git `

   navigate to the downloaded directory:

   ` cd Mshdist `

   then create build directory and compile the project using cmake
   ```
   mkdir build
   cd build
   cmake ..
   make
   make install
   ```

#### Usage

* A short documentation is supplied (see the [documentation](documentation) folder).

#### Quickstart

*

#### Authors & contributors
* mshdist is developped and maintained by Charles Dapogny (Université Joseph Fourier) and Pascal Frey (Université Pierre et Marie Curie).
* Contributors to this project are warmly welcome.

#### License
* mshdist is given under the [terms of the GNU Lesser General Public License] (LICENSE.md).

* If you use mshdist in your work, please refer to the journal article:

C. Dapogny, P. Frey, _Computation of the signed distance function to a discrete contour on adapted triangulation_, Calcolo, Volume 49, Issue 3, pp.193-219 (2012).
