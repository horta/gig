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
