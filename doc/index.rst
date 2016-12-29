===================
Gig's documentation
===================


Draws from the Generalized Inverse Gaussian distribution

.. math::

  f(x) = x^{\lambda - 1} e^{-\frac{1}{2}(\chi/x + \psi x)}

where :math:`\lambda` is any real number, :math:`\chi` must be nonnegative
(nonpositive) for positive (negative) :math:`\lambda` and :math:`\psi` must be
be nonnegative (nonpositive) for negative (positive) :math:`\lambda`.

-------
Install
-------

In the project folder, type

.. code-block:: bash

  mkdir build
  cd build
  cmake ..
  make

It should create the shared and static libraries

.. code-block:: bash

  libgig.[version].[extension]
  libgig_static.[extension]

You can enter

.. code-block:: bash

  make test

to test it and

.. code-block:: bash

  make install

to install.

-------------
Usage example
-------------

Suppose you have the file

.. code-block:: cpp

  /* example.cpp */
  #include "gig/gig.h"

  #include <iostream>

  int main()
  {
    Random random(1);

    double lambda = 2.1;
    double chi = 0.1;
    double psi = 1.0;

    std::cout << random.gig(lambda, chi, psi) << std::endl;
  }

Compiling, linking, and running it via

.. code-block:: bash

  cpp libgig_static.a example.cpp -o example
  ./example

should print::

  -2.049961
