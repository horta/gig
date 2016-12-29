===================
Gig's documentation
===================


Draws from the Generalized Inverse Gaussian distribution

.. math::

  f(x) = x^{\lambda - 1} e^{-\frac{1}{2}(\chi/x + \psi x)}

where :math:`\lambda` is any real number, :math:`\chi` must be be nonnegative
for positive :math:`\lambda` or nonpositive otherwise and :math:`\psi` must be
be nonnegative for negative :math:`\lambda` or nonnegative otherwise.

--------
Build it
--------

In the project folder, type

.. code-block:: bash

  mkdir build
  cd build
  cmake ..
  make

It should create the library

.. code-block:: bash

  gig/libgig.[extension]

And if you want to test it:

.. code-block:: bash

  make test
