# Mshdist

| Branch      |  Matrix  |
|:-----------:|:--------:|
| master      | [![Build](https://github.com/cbritopacheco/Mshdist/actions/workflows/Build.yml/badge.svg?branch=master)](https://github.com/cbritopacheco/Mshdist/actions/workflows/Build.yml) |
| develop     | [![Build](https://github.com/cbritopacheco/Mshdist/actions/workflows/Build.yml/badge.svg?branch=develop)](https://github.com/cbritopacheco/Mshdist/actions/workflows/Build.yml) |

## About
Mshdist is a simple algorithm to generate the signed distance function to given
objects in two and three space dimensions.

This repo is a fork of the [original
Mshdist](https://github.com/ISCDtoolbox/Mshdist) to be used as a dependency in
[Rodin](https://github.com/cbritopacheco/rodin).

## Building

```bash
git clone --recursive https://github.com/cbritopacheco/Mshdist
cd Mshdist
mkdir build && cd build
cmake ..
make -j4
```

## License

mshdist is given under the [terms of the GNU Lesser General Public License](LICENSE.md).

## Credits
Mshdist is originally based on the paper: [Dapogny, Charles, and Pascal Frey. "Computation of the signed distance function to a discrete contour on adapted triangulation." Calcolo 49.3 (2012): 193-219.](https://link.springer.com/article/10.1007/s10092-011-0051-z)
