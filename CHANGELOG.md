# Changelog

changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## Unreleased

### Added

 - Use new Info system of PaRSEC to register cublas and
   CuDnSolver handles and allow to run all kernels of
   POTRF on GPUs
 - zgemm uses improved control flow to keep active set in GPU memory when
   the problem size does not hold in GPU memory (see Scala'19).
 - Support for potri functions (trtri+lauum), and corresponding testings.

### Changed

 - Improved ztrmm with all the matrix reads unified.

### Deprecated

### Removed

 - DPLASMA does not depend anymore on external COREBLAS or PLASMA, but only on BLAS
   and LAPACKE

### Fixed

  - Fix bug in symmetric/hermitian norms when `A.mt < P`.

### Security


## Prior versions

This list contains the main features as well as overviews of specific
bug fixes (and other actions) for each version of DPLASMA since
inception (basically the split from the PaRSEC project).

 - Split DPLASMA and PaRSEC into separate repositories. PaRSEC moves from
   cmake-2.0 to cmake-3.12, using targets. Targets are exported for
   third-party integration
 - Add the doxygen documentation generation to the configure system.
 - Add support in the runtime for user-defined properties evaluated at
   runtime and easy to export through a shared memory region (see: PR
   229 visualization-tools)
 - Change several default parameters for the DPLASMA tests.
