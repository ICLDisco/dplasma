DPLASMA_ is the leading implementation of a dense linear algebra package
for distributed, accelerated, heterogeneous systems. It is designed to 
deliver sustained performance for distributed systems where each node featuring
multiple sockets of multicore processors, and if available, accelerators like
GPUs or Intel Xeon Phi. DPLASMA achieves this objective through the state of
the art PaRSEC_ runtime, porting the Parallel Linear Algebra Software for
Multicore Architectures (PLASMA) algorithms to the distributed memory realm.

.. _DPLASMA: https://github.com/icldisco/dplasma
.. _PaRSEC: https://github.com/icldisco/parsec

.. contents:: Table of Contents


Features
========

DPLASMA’s feature set includes:

* Linear systems of equations (Cholesky, LU [inc. pivoting, PP], LDL [prototype]).
* Least squares (QR & LQ).
* Symmetric eigenvalue problem (Reduction to Band [prototype]).
* Level 3 Tile BLAS (GEMM, TRSM, TRMM, HEMM/SYMM, HERK/SYRK, HER2K/SYR2k.
* Covers double real, double complex, single real, and single complex (D, Z, S, C) precisions.
* Provides ScaLAPACK-compatible interface for matrices in F77 column-major layout
* Supports: Linux, Windows, Mac OS X, UN*X (depends on MPI, hwloc)


History
=======

DPLASMA started as the distributed version of ICL PLASMA linear algebra package
implemented on top of the PaRSEC_ runtime, a generic framework for architecture
aware scheduling and management of micro-tasks on distributed many-core heterogeneous
architectures. For a long period, the DPLASMA library was packaged as part of the PaRSEC
project. DPLASMA is now an independent project and uses PaRSEC as a dependency.


Installation
============

Please read the INSTALL.rst_ file for the software dependencies and the installation
instructions. In addition, the `PaRSEC wiki`_ might contain additional information.

.. _INSTALL.rst: https://github.com/ICLDisco/dplasma/blob/master/INSTALL.rst
.. _`PaRSEC wiki`: https://github.com/icldisco/parsec/wiki


Links
=====

When referencing DPLASMA in a publication please use the following reference:
`[10.1109/IPDPS.2011.299]`_: Flexible Development of Dense Linear Algebra Algorithms on Massively Parallel Architectures with DPLASMA [bibtex_].

.. _[10.1109/IPDPS.2011.299]: http://www.icl.utk.edu/node/617
.. _bibtex: http://www.icl.utk.edu/publications/export/bibtex/617
 
If you are interested in this project and want to join the users community, our
mailman will be happy to accept you on the project user (moderated) mailing list
at dplasma-users@icl.utk.edu.

Happy hacking,
  The PaRSEC team.

