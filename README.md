# Mshdist

Mshdist is a simple algorithm to generate the signed distance function to given
objects in two and three space dimensions.

This repo is a fork of the [original
Mshdist](https://github.com/ISCDtoolbox/Mshdist) to be used as a dependency in
[Rodin](https://github.com/cbritopacheco/rodin).

# Building

```bash
git clone --recursive https://github.com/cbritopacheco/Mshdist
cd Mshdist
mkdir build && cd build
cmake ..
make -j4
```

# License

mshdist is given under the [terms of the GNU Lesser General Public License] (LICENSE.md).

If you use mshdist in your work, please refer to the journal article:

C. Dapogny, P. Frey, _Computation of the signed distance function to a discrete
contour on adapted triangulation_, Calcolo, Volume 49, Issue 3, pp.193-219
(2012).
