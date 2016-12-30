# gig

[![Travis](https://img.shields.io/travis/Horta/gig.svg?style=flat-square)](https://travis-ci.org/Horta/gig)
[![Documentation Status](https://readthedocs.org/projects/gig/badge/?style=flat-square&version=latest)](http://gig.readthedocs.io/en/latest/?badge=latest)

Draw sample from the Generalized Inverse Gaussian distribution.

## Install

Download and unpack the [latest version](https://github.com/Horta/gig/releases/latest).
In the unpacked folder, type

```bash
mkdir build
cd build
cmake ..
make
```

It should create the shared and static libraries

```bash
libgig.[version].[extension]
libgig_static.[extension]
```

You can enter

```bash
make test
```

to test it and

```bash
make install
```

to install.

## Usage example

Suppose you have the file

```cpp
/* example.cpp */
#include "gig/gig.h"

#include <random>
#include <iostream>

int main()
{
  Random<std::default_random_engine> random(1);

  double lambda = 2.1;
  double chi = 0.1;
  double psi = 1.0;

  std::cout << random.gig(lambda, chi, psi) << std::endl;
}
```

Compiling, linking, and running it via

```bash
g++ libgig_static.a  example.cpp -o example
./example
```

should print

```
1.30869
```

## Disclaimer

This library is simply a wrapper around Josef Leydold and Wolfgang Hormann's
implementation of a GIG sampler found in the
[GIGrvg package](https://cran.r-project.org/web/packages/GIGrvg/GIGrvg.pdf).

## License

This project is licensed under the MIT License - see the
[LICENSE](LICENSE) file for details
