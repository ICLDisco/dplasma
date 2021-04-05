==================
INSTALLING DPLASMA
==================

.. contents:: Table of Contents

Software Dependencies
=====================

To compile DPLASMA on a new platform, you will need some of the software
below. From 1 to 4 (included) they are mandatory. Everything else is
optional, they provide nice features not critical to the normal usage
of this software package.

1. cmake version 3.12 or above. cmake can be found in the debian
   package cmake, or as sources at the CMake_ download page
2. A version of the Basic Linear Algebra Subroutines (BLAS) with LAPACKE_
   (or an external LAPACK/LAPACKE to completement your BLAS, see below)
3. Any MPI library Open MPI, MPICH2, MVAPICH or any vendor blessed
   implementation.
4. The PaRSEC_ Runtime
5. Optional (but highly recommended): hwloc_ for processor and memory 
   locality features in PaRSEC
6. Optional: CUDA_ for Nvidia hardware acceleration

.. _CMake: http://www.cmake.org/
.. _LAPACKE: https://github.com/Reference-LAPACK/lapack
.. _PaRSEC: https://bitbucket.org/icldistcomp/parsec/
.. _hwloc: http://www.open-mpi.org/projects/hwloc/
.. _CUDA: https://developer.nvidia.com/cuda-zone


Configuring DPLASMA for a new platform
======================================

DPLASMA is a CMake_ built project. CMake has a comparable goal to
configure_, but it's subtly different. For one thing, CMake display the
commands with colors, but this is not necessarily its most prominent
feature.

CMake keeps everything it found hitherto in a cache file named
``CMakeCache.txt``. Until you have successfully configured DPLASMA,
remove the ``CMakeCache.txt`` file each time you run ``cmake``.

.. _configure: https://www.gnu.org/software/autoconf/

The configure script
--------------------

Passing options to CMake can be confusing. For that reason we have
designed the ``configure`` script that basically wraps around the
invocation of ``cmake`` with a tested and trusted feel to it.

.. code:: bash

  configure --prefix=$INSTALLDIR --with-mpi=$MPI_DIR --with-cuda --disable-debug

will produce a ``cmake`` command line that matches these options,
and execute it.

You can review what ``cmake`` invocation has been produced by looking
into ``config.log``.

It also produces a ``config.status`` script that helps redo the last
``configure`` step, while ``config.status --recheck`` will also clean
the CMake caches.

Not all options you can pass to PaRSEC exist as a ``--enable-xxx``
``--with-yyy`` configure argument. You can pass environment variables
to the produced ``cmake`` command as well as CMake *defines* (both
will appear in the ``config.log``) by using the following form:

.. code:: bash

    configure CC=icc FC=ftn CXX=icpc -DPARSEC_EAGER_LIMIT=0

Plarform Files
--------------

Platform files, found in ``contrib/platforms`` let us distribute recipes
for well known systems that may be similar to a supercomputer near you.
For example, the ``ibm.ac922.summit`` file is intended to compile on the
eponyme Oak Ridge Leadership Computing Facility system.

.. code:: bash

  configure --prefix=$INSTALL_DIR --with-platform=ibm.ac922.summit --disable-debug

This call should get you running in no time on that machine, and you
may still customize and overide the platform file with command line
arguments.

We also provide a ``macosx`` platform file that helps dealing with the
detection of the Fortran compiler on this architecture.

Of course you may edit and produce your own platform files for your
favorite computer. These are shell script that execute in the context
of the main configure script. For example, our continuous integration
system is named *saturn*, in that script you will find examples of
how one sets some default options.

.. code:: bash

  with_hwloc=${HWLOC_ROOT:="/spack/opt/spack/linux-scientific7-x86_64/gcc-7.3.0/hwloc-1.11.11-nu65xwuyodswr74llx3ymi67hgd6vmwe"}

  # BLAS: use MKL
  [ -z "${MKLROOT}" ] || module load intel-mkl/2019.3.199/gcc-7.3.0-2pn4
  with_blas=Intel10_64lp_seq

  # Slurm test options
  CMAKE_DEFINES+=" -DCTEST_MPI_LAUNCHER=\"srun -Ccauchy -N\" -DCTEST_SHM_LAUNCHER=\"srun -Ccauchy\" -DCTEST_GPU_LAUNCHER_OPTIONS=-Cgtx1060"

As you can see, the platform file may contain commands, shell scripts,
load environment modules_, etc. Of note are the ``CMAKE_DEFINES`` and
``ENVVARS`` variables which control what ``-DX=Y`` options are appended
, and ``A=B`` environment are prepended to the ``cmake`` invocation,
respectively.

Submodule or External PaRSEC
----------------------------

By default, DPLASMA will try to detect as system (or speficied in the
``PaRSEC_ROOT`` environment variable) automatically. If an installed
PaRSEC is not found, DPLASMA will download an appropriate version of
PaRSEC from ``bitbucket.org`` and setup a ``git submodule``. This
Submodule PaRSEC will be configured and built at the same time as
DPLASMA. Passing ``--without-parsec`` to ``configure``  will force using
the submodule PaRSEC instead of looking for an installed version.

Conversely, you can prevent loading the Submodule PaRSEC by setting
``--with-parsec``. You can select a particular externally installed
PaRSEC by setting the configure option 
``--with-parsec=$PARSEC_INSTALL_DIRECTORY``.

Note that many of the ``configure`` options apply only to the submodule
PaRSEC and have no effect when you are using an external PaRSEC. Setting
these will result in a warning by CMake that some variables have been 
defined but unused.

Cross Compiling
---------------

On some system, the build machine cannot execute the code produced for
compute nodes. An example is the ANL Theta system, a Cray XC40
with Xeon Phi nodes and Haswell build frontends.

Cross compiling is heavily reliant on the *platform file* feature.
For example, on the Theta system, one can cross compile by simply
calling

.. code:: bash

  configure --with-platform=cray.xc40.theta

In this case, the configuration stage will also include a build stage
to produce some of the utilities needed to compile PaRSEC. After
the configure state has completed, you will find in your build directory
a subdirectory named ``native`` that contains profiling and devellopper
tools that can be used on the frontend system.

After the configure step has completed, the build step is carried out
as usual by simply using ``make``.

If you face a new system where you need to cross compile, a good start
is to copy the ``contrib/platforms/cray.xc40.theta`` file, and
customize it according to your needs.

Note that you will most probably need to produce your own ``toolchain``
CMake cross-compilation file. More information can be found about them
on the cmake-toolchain_ web page.

.. _cmake-toolchain: https://cmake.org/cmake/help/v3.14/manual/cmake-toolchains.7.html?highlight=cross

Legacy Configurations
---------------------

Of course, you can always directly invoke ``cmake``. You can take
inspiration from the command produced from the ``configure`` script,
or you can look at the obsolete ``contrib/platforms/legacy/config.inc``.

.. code:: bash

  rm -f CMakeCache.txt
  cmake . -G 'Unix Makefiles' -DPARSEC_DIST_WITH_MPI=ON

``contrib/platforms/legacy`` also contains shell scripts that we used to
configure on older systems. ``config.jaguar`` is for, you got it, XT5,
etc. If your system is similar to one of these old systems, we advise
you to start from a modern platform file and tweak from there by importing
the content of the old scripts. Unlike modern platform files, legacy
scripts are shell scripts that can be executed directly from desired
build directory (VPATH or not).


Full configuration example
--------------------------

Hopefully, once the expected arguments are provided the output will look similar to

.. code:: console

  ### This program was invoked with the following command line
  #
      ../dplasma/configure  --with-platform=ibm.ac922.summit --enable-debug=noisier\,paranoid
  #
  #################################################
  # Platform ibm.ac922.summit
  # This file is for a compilation on OLCF Summit.
  #   configure --with-platform=ibm.ac922.summit ...
  # Set preferences and dependencies for the
  # ibm.ac922.summit system executables and libs
  #   CC=mpicc CXX=mpiCC FC=mpif90
  #
  
  The following have been reloaded with a version change:
    1) cmake/3.14.2 => cmake/3.15.2
  
  ### CMake generated invocation
  #
       LAPACKE_ROOT=/ccs/home/bouteilla/parsec/dplasma/lapack CC=mpicc CXX=mpicxx FC=mpif90 CFLAGS='' LDFLAGS='' /autofs/nccs-svm1_sw/summit/.swci/0-core/opt/spack/20180914/linux-rhel7-ppc64le/gcc-4.8.5/cmake-3.15.2-xit2o3iepxvqbyku77lwcugufilztu7t/bin/cmake -G 'Unix Makefiles' /ccs/home/bouteilla/parsec/summit.debug.dplasma/../dplasma  -DBLAS_LIBRARIES='/sw/summit/essl/6.2.0-20190419/essl/6.2/lib64/libessl.so' -DBLA_VENDOR=IBMESSL -DCMAKE_INSTALL_PREFIX=/usr/local -DCMAKE_BUILD_TYPE=Debug -DPARSEC_DEBUG_PARANOID=ON -DPARSEC_DEBUG_NOISIER=ON -DPARSEC_GPU_WITH_CUDA=ON
  #
  Removing Cmake Cache...
  -- The C compiler identification is XLClang 16.1.1.3
  -- Check for working C compiler: /autofs/nccs-svm1_sw/summit/.swci/1-compute/opt/spack/20180914/linux-rhel7-ppc64le/xl-16.1.1-3/spectrum-mpi-10.3.0.1-20190611-aqjt3jo53mogrrhcrd2iufr435azcaha/bin/mpicc
  -- Check for working C compiler: /autofs/nccs-svm1_sw/summit/.swci/1-compute/opt/spack/20180914/linux-rhel7-ppc64le/xl-16.1.1-3/spectrum-mpi-10.3.0.1-20190611-aqjt3jo53mogrrhcrd2iufr435azcaha/bin/mpicc -- works
  -- Detecting C compiler ABI info
  -- Detecting C compiler ABI info - done
  -- Detecting C compile features
  -- Detecting C compile features - done
  -- The Fortran compiler identification is XL 16.1.1
  -- Check for working Fortran compiler: /autofs/nccs-svm1_sw/summit/.swci/1-compute/opt/spack/20180914/linux-rhel7-ppc64le/xl-16.1.1-3/spectrum-mpi-10.3.0.1-20190611-aqjt3jo53mogrrhcrd2iufr435azcaha/bin/mpif90
  -- Check for working Fortran compiler: /autofs/nccs-svm1_sw/summit/.swci/1-compute/opt/spack/20180914/linux-rhel7-ppc64le/xl-16.1.1-3/spectrum-mpi-10.3.0.1-20190611-aqjt3jo53mogrrhcrd2iufr435azcaha/bin/mpif90  -- works
  -- Detecting Fortran compiler ABI info
  -- Detecting Fortran compiler ABI info - done
  -- Checking whether /autofs/nccs-svm1_sw/summit/.swci/1-compute/opt/spack/20180914/linux-rhel7-ppc64le/xl-16.1.1-3/spectrum-mpi-10.3.0.1-20190611-aqjt3jo53mogrrhcrd2iufr435azcaha/bin/mpif90 supports Fortran 90
  -- Checking whether /autofs/nccs-svm1_sw/summit/.swci/1-compute/opt/spack/20180914/linux-rhel7-ppc64le/xl-16.1.1-3/spectrum-mpi-10.3.0.1-20190611-aqjt3jo53mogrrhcrd2iufr435azcaha/bin/mpif90 supports Fortran 90 -- yes
  -- Detecting Fortran/C Interface
  -- Detecting Fortran/C Interface - Found GLOBAL and MODULE mangling
  -- Found BLAS: /sw/summit/essl/6.2.0-20190419/essl/6.2/lib64/libessl.so
  -- Looking for zgemm
  -- Looking for zgemm - found
  -- Looking for Fortran zgeqrf
  -- Looking for Fortran zgeqrf - found
  -- Performing Test BLAS_HAS_CBLAS
  -- Performing Test BLAS_HAS_CBLAS - Success
  -- Performing Test BLAS_HAS_LAPACKE
  -- Performing Test BLAS_HAS_LAPACKE - Success
  -- Found LAPACKE: /sw/summit/essl/6.2.0-20190419/essl/6.2/lib64/libessl.so  found components:  BLAS CBLAS LAPACK LAPACKE
  -- Found LAPACKE and defined the following imported targets:
  --   - LAPACKE::LAPACKE:
  --       + include:      /sw/summit/essl/6.2.0-20190419/essl/6.2/include;/ccs/home/bouteilla/parsec/dplasma/lapack/LAPACKE/include
  --       + library:      /sw/summit/essl/6.2.0-20190419/essl/6.2/lib64/libessl.so
  --       + dependencies: /ccs/home/bouteilla/parsec/dplasma/lapack/liblapacke.a;
  --   - LAPACKE::LAPACK:
  --       + include:      /sw/summit/essl/6.2.0-20190419/essl/6.2/include;/ccs/home/bouteilla/parsec/dplasma/lapack/LAPACKE/include
  --       + library:      /sw/summit/essl/6.2.0-20190419/essl/6.2/lib64/libessl.so
  --       + dependencies: /ccs/home/bouteilla/parsec/dplasma/lapack/liblapack.a;
  --   - LAPACKE::CBLAS:
  --       + include:      /sw/summit/essl/6.2.0-20190419/essl/6.2/include;/ccs/home/bouteilla/parsec/dplasma/lapack/LAPACKE/include
  --       + library:      /sw/summit/essl/6.2.0-20190419/essl/6.2/lib64/libessl.so
  --       + dependencies:
  --   - LAPACKE::BLAS:
  --       + include:      /sw/summit/essl/6.2.0-20190419/essl/6.2/include;/ccs/home/bouteilla/parsec/dplasma/lapack/LAPACKE/include
  --       + library:      /sw/summit/essl/6.2.0-20190419/essl/6.2/lib64/libessl.so
  --       + dependencies: /sw/summit/essl/6.2.0-20190419/essl/6.2/lib64/libessl.so;
  -- Looking for timersub
  -- Looking for timersub - found
  -- Looking for asprintf
  -- Looking for asprintf - not found
  -- Looking for asprintf
  -- Looking for asprintf - found
  -- Found PythonInterp: /usr/bin/python (found version "2.7.5")
  -- ########################################################################
  -- #             Configuring internal submodule PaRSEC runtime!
  -- The CXX compiler identification is XLClang 16.1.1.3
  -- Check for working CXX compiler: /autofs/nccs-svm1_sw/summit/.swci/1-compute/opt/spack/20180914/linux-rhel7-ppc64le/xl-16.1.1-3/spectrum-mpi-10.3.0.1-20190611-aqjt3jo53mogrrhcrd2iufr435azcaha/bin/mpicxx
  -- Check for working CXX compiler: /autofs/nccs-svm1_sw/summit/.swci/1-compute/opt/spack/20180914/linux-rhel7-ppc64le/xl-16.1.1-3/spectrum-mpi-10.3.0.1-20190611-aqjt3jo53mogrrhcrd2iufr435azcaha/bin/mpicxx -- works
  -- Detecting CXX compiler ABI info
  -- Detecting CXX compiler ABI info - done
  -- Detecting CXX compile features
  -- Detecting CXX compile features - done
  -- Found BISON: /usr/bin/bison (found version "3.0.4")
  -- Found FLEX: /usr/bin/flex (found version "2.5.37")
  -- Building for target ppc64le
  -- Found target for PPC
  -- Performing Test C_M32or64
  -- Performing Test C_M32or64 - Success
  -- Performing Test PARSEC_HAVE_STD_C1x
  -- Performing Test PARSEC_HAVE_STD_C1x - Success
  -- Performing Test PARSEC_HAVE_STD_C99
  -- Performing Test PARSEC_HAVE_STD_C99 - Success
  -- Performing Test PARSEC_HAVE_WD
  -- Performing Test PARSEC_HAVE_WD - Failed
  -- Performing Test PARSEC_HAVE_G3
  -- Performing Test PARSEC_HAVE_G3 - Success
  -- Looking for sys/types.h
  -- Looking for sys/types.h - found
  -- Looking for stdint.h
  -- Looking for stdint.h - found
  -- Looking for stddef.h
  -- Looking for stddef.h - found
  -- Check size of __int128_t
  -- Check size of __int128_t - done
  -- Performing Test PARSEC_COMPILER_C11_COMPLIANT
  -- Performing Test PARSEC_COMPILER_C11_COMPLIANT - Failed
  -- Performing Test PARSEC_ATOMIC_USE_GCC_32_BUILTINS
  -- Performing Test PARSEC_ATOMIC_USE_GCC_32_BUILTINS - Success
  -- Performing Test PARSEC_ATOMIC_USE_GCC_64_BUILTINS
  -- Performing Test PARSEC_ATOMIC_USE_GCC_64_BUILTINS - Success
  -- Performing Test PARSEC_ATOMIC_USE_GCC_128_BUILTINS
  -- Performing Test PARSEC_ATOMIC_USE_GCC_128_BUILTINS - Failed
  -- Performing Test PARSEC_ATOMIC_USE_GCC_128_BUILTINS
  -- Performing Test PARSEC_ATOMIC_USE_GCC_128_BUILTINS - Failed
  -- Performing Test PARSEC_ATOMIC_USE_XLC_32_BUILTINS
  -- Performing Test PARSEC_ATOMIC_USE_XLC_32_BUILTINS - Success
  -- Performing Test PARSEC_ATOMIC_USE_XLC_64_BUILTINS
  -- Performing Test PARSEC_ATOMIC_USE_XLC_64_BUILTINS - Success
  -- Performing Test PARSEC_ATOMIC_USE_XLC_LLSC_32_BUILTINS
  -- Performing Test PARSEC_ATOMIC_USE_XLC_LLSC_32_BUILTINS - Success
  -- Performing Test PARSEC_ATOMIC_USE_XLC_LLSC_64_BUILTINS
  -- Performing Test PARSEC_ATOMIC_USE_XLC_LLSC_64_BUILTINS - Success
  -- Performing Test PARSEC_ATOMIC_USE_MIPOSPRO_32_BUILTINS
  -- Performing Test PARSEC_ATOMIC_USE_MIPOSPRO_32_BUILTINS - Failed
  -- Performing Test PARSEC_ATOMIC_USE_SUN_32
  -- Performing Test PARSEC_ATOMIC_USE_SUN_32 - Failed
  --       support for 32 bits atomics - found
  --       support for 64 bits atomics - found
  --       support for XL LL/SC atomics - found
  -- Looking for pthread.h
  -- Looking for pthread.h - found
  -- Performing Test CMAKE_HAVE_LIBC_PTHREAD
  -- Performing Test CMAKE_HAVE_LIBC_PTHREAD - Success
  -- Found Threads: TRUE
  -- Looking for pthread_getspecific
  -- Looking for pthread_getspecific - found
  -- Looking for pthread_barrier_init
  -- Looking for pthread_barrier_init - found
  -- Looking for sched_setaffinity
  -- Looking for sched_setaffinity - found
  -- Performing Test PARSEC_HAVE_TIMESPEC_TV_NSEC
  -- Performing Test PARSEC_HAVE_TIMESPEC_TV_NSEC - Success
  -- Looking for clock_gettime in c
  -- Looking for clock_gettime in c - found
  -- Looking for include file stdarg.h
  -- Looking for include file stdarg.h - found
  -- Performing Test PARSEC_HAVE_VA_COPY
  -- Performing Test PARSEC_HAVE_VA_COPY - Success
  -- Performing Test PARSEC_HAVE_ATTRIBUTE_FORMAT_PRINTF
  -- Performing Test PARSEC_HAVE_ATTRIBUTE_FORMAT_PRINTF - Success
  -- Performing Test PARSEC_HAVE_THREAD_LOCAL
  -- Performing Test PARSEC_HAVE_THREAD_LOCAL - Success
  -- Looking for asprintf
  -- Looking for asprintf - found
  -- Looking for vasprintf
  -- Looking for vasprintf - found
  -- Looking for include file unistd.h
  -- Looking for include file unistd.h - found
  -- Looking for include file getopt.h
  -- Looking for include file getopt.h - found
  -- Looking for getopt_long
  -- Looking for getopt_long - found
  -- Looking for include file errno.h
  -- Looking for include file errno.h - found
  -- Looking for include file stddef.h
  -- Looking for include file stddef.h - found
  -- Looking for include file stdbool.h
  -- Looking for include file stdbool.h - found
  -- Looking for include file ctype.h
  -- Looking for include file ctype.h - found
  -- Performing Test PARSEC_HAVE_BUILTIN_CPU
  -- Performing Test PARSEC_HAVE_BUILTIN_CPU - Failed
  -- Looking for getrusage
  -- Looking for getrusage - found
  -- Looking for RUSAGE_THREAD
  -- Looking for RUSAGE_THREAD - not found
  -- Looking for RUSAGE_THREAD
  -- Looking for RUSAGE_THREAD - found
  -- Looking for include file limits.h
  -- Looking for include file limits.h - found
  -- Looking for include file string.h
  -- Looking for include file string.h - found
  -- Looking for include file libgen.h
  -- Looking for include file libgen.h - found
  -- Looking for include file complex.h
  -- Looking for include file complex.h - found
  -- Looking for include file sys/param.h
  -- Looking for include file sys/param.h - found
  -- Looking for include file sys/types.h
  -- Looking for include file sys/types.h - found
  -- Looking for include file syslog.h
  -- Looking for include file syslog.h - found
  -- Performing Test PARSEC_HAVE_ATTRIBUTE_ALWAYS_INLINE
  -- Performing Test PARSEC_HAVE_ATTRIBUTE_ALWAYS_INLINE - Success
  -- Performing Test PARSEC_HAVE_ATTRIBUTE_VISIBILITY
  -- Performing Test PARSEC_HAVE_ATTRIBUTE_VISIBILITY - Success
  -- Performing Test PARSEC_HAVE_BUILTIN_EXPECT
  -- Performing Test PARSEC_HAVE_BUILTIN_EXPECT - Success
  -- Found HWLOC: /usr/lib64/libhwloc.so
  -- Performing Test PARSEC_HAVE_HWLOC_PARENT_MEMBER
  -- Performing Test PARSEC_HAVE_HWLOC_PARENT_MEMBER - Success
  -- Performing Test PARSEC_HAVE_HWLOC_CACHE_ATTR
  -- Performing Test PARSEC_HAVE_HWLOC_CACHE_ATTR - Success
  -- Performing Test PARSEC_HAVE_HWLOC_OBJ_PU
  -- Performing Test PARSEC_HAVE_HWLOC_OBJ_PU - Success
  -- Looking for hwloc_bitmap_free in /usr/lib64/libhwloc.so
  -- Looking for hwloc_bitmap_free in /usr/lib64/libhwloc.so - found
  -- Performing Test CC_CONTAINS_MPI
  -- Performing Test CC_CONTAINS_MPI - Success
  -- Looking for MPI_Type_create_resized
  -- Looking for MPI_Type_create_resized - found
  -- Performing Test PARSEC_HAVE_MPI_OVERTAKE
  -- Performing Test PARSEC_HAVE_MPI_OVERTAKE - Success
  -- Found CUDA: /sw/summit/cuda/10.1.168 (found version "10.1")
  -- Found CUDA 10.1 in /sw/summit/cuda/10.1.168
  -- Looking for cudaDeviceCanAccessPeer
  -- Looking for cudaDeviceCanAccessPeer - found
  -- Add -q64 and -nofor_main to the Fortran linker.
  CMAKE_Fortran_COMPILER full path: /autofs/nccs-svm1_sw/summit/.swci/1-compute/opt/spack/20180914/linux-rhel7-ppc64le/xl-16.1.1-3/spectrum-mpi-10.3.0.1-20190611-aqjt3jo53mogrrhcrd2iufr435azcaha/bin/mpif90
  Fortran compiler: mpif90
  No optimized Fortran compiler flags are known, we just try -O2...
  -- Checking for module 'libgvc'
  --   No package 'libgvc' found
  -- Could NOT find GRAPHVIZ (missing: GRAPHVIZ_LIBRARY GRAPHVIZ_INCLUDE_DIR)
  -- Could NOT find Cython (missing: CYTHON_EXECUTABLE) (Required is at least version "0.21.2")
  -- Looking for shm_open
  -- Looking for shm_open - not found
  -- Looking for shm_open in rt
  -- Looking for shm_open in rt - found
  -- PARSEC Modular Component Architecture (MCA) discovery:
  -- -- Found Component `pins'
  -- Module alperf not selectable: PARSEC_PROF_TRACE disabled.
  -- ---- Module `iterators_checker' is ON
  -- Module papi not selectable: PARSEC_PROF_TRACE disabled.
  -- ---- Module `print_steals' is ON
  -- ---- Module `ptg_to_dtd' is ON
  -- Module task_profiler not selectable: PARSEC_PROF_TRACE disabled.
  -- Component pins sources: mca/pins/pins.c;mca/pins/pins_init.c
  -- -- Found Component `sched'
  -- ---- Module `ap' is ON
  -- ---- Module `gd' is ON
  -- ---- Module `ip' is ON
  -- ---- Module `lfq' is ON
  -- ---- Module `lhq' is ON
  -- ---- Module `ll' is ON
  -- ---- Module `ltq' is ON
  -- ---- Module `pbq' is ON
  -- ---- Module `rnd' is ON
  -- ---- Module `spq' is ON
  -- Component sched sources:
  -- PARSEC Modular Component Architecture (MCA) discovery done.
  -- Looking for PARSEC_ATOMIC_HAS_ATOMIC_CAS_INT128
  -- Looking for PARSEC_ATOMIC_HAS_ATOMIC_CAS_INT128 - not found
  -- Check size of ((parsec_lifo_t*)0)->lifo_head
  -- Check size of ((parsec_lifo_t*)0)->lifo_head - done
  -- Internal PaRSEC does not use CAS on int128_t. Keeping parsec_options.h unchanged
  
  
  Configuration flags:
    CMAKE_C_FLAGS          =  -q64 -qlanglvl=extc99
    CMAKE_C_LDFLAGS        =  -q64
    CMAKE_EXE_LINKER_FLAGS =
    EXTRA_LIBS             = /usr/lib64/libhwloc.so
  
  
  
  -- #             Configuring internal submodule PaRSEC runtime: DONE!
  -- ########################################################################
  -- CUDA support for DPLASMA enabled
  -- Looking for include file complex.h
  -- Looking for include file complex.h - found
  -- Generate precision dependencies in /ccs/home/bouteilla/parsec/dplasma/include        generated_headers
  -- Generate precision dependencies in /ccs/home/bouteilla/parsec/dplasma/cores  generated_headers
  -- Generate precision dependencies in /ccs/home/bouteilla/parsec/dplasma/cores  generated_files
  -- Generate precision dependencies in /ccs/home/bouteilla/parsec/dplasma/cores  all_precisions_files
  -- Generate precision dependencies in /ccs/home/bouteilla/parsec/dplasma/cores  cplx_files
  -- Generate precision dependencies in /ccs/home/bouteilla/parsec/dplasma/cores  generated_cuda_files
  -- Generate precision dependencies in /ccs/home/bouteilla/parsec/dplasma/lib    generated_jdf
  -- Generate precision dependencies in /ccs/home/bouteilla/parsec/dplasma/lib    generated_wrappers
  -- Generate precision dependencies in /ccs/home/bouteilla/parsec/dplasma/tests  generated_testings
  -- Configuring done
  -- Generating done
  -- Build files have been written to: /ccs/home/bouteilla/parsec/summit.debug.dplasma

If this is done, congratulations, DPLASMA is configured and you're ready for
building and testing the system.

Troubleshooting
---------------

In the unlikely case something goes wrong, read carefully the error message. We
spend a significant amount of time trying to output something meaningful for you
and for us (in case you need help to debug/understand). If the output is not
helpful enough to fix the problem, you should contact us via the PaRSEC user
mailing list and provide the CMake command and the flags, the output as well as
the files CMakeFiles/CMakeError.log and CMakeFiles/CMakeOutput.log.

We use quite a few packages that are optional, don't panic if they are not found
during the configuration. However, some of them are critical for increasing the
performance (such as HWLOC).

Check that you have a working MPI somewhere accessible (``mpicc`` and ``mpirun`` should
be in your PATH, except on Cray systems where you should use the ``cc`` wrapper).

If you have strange behavior, check that you have a success line for one of the
possible atomic backends that make sense for your local environment (i.e.,
C11 or GNU atomics depending on GCC versions, XLC on BlueGene machines, etc.).
If not, the atomic operations will not work, and that is damageable for the good
operation of PaRSEC. Note how in the shown configuration below, it takes
several attempts to get the right flags to use 128 bits atomic operations, but
in the end all looks good here.

.. code:: console

  -- Found target X86_64
  ...
  -- Performing Test PARSEC_ATOMIC_USE_C11_128
  -- Performing Test PARSEC_ATOMIC_USE_C11_128 - Failed
  -- Performing Test PARSEC_ATOMIC_USE_C11_128
  -- Performing Test PARSEC_ATOMIC_USE_C11_128 - Failed
  -- Performing Test PARSEC_ATOMIC_USE_C11_128
  -- Performing Test PARSEC_ATOMIC_USE_C11_128 - Success
  --       support for 32 bits atomics - found
  --       support for 64 bits atomics - found
  --       support for 128 bits atomics - found

CMake behavior can be modified from what your environment variables contain.
For example environment modules_, a popular way to load software on Cray,
DOE and NERSC supercomputers, can set many variables that will change the
outcome of the CMake configuration stage.

CC
  to choose your C compiler
CFLAGS
  to change your C compilation flags
LDFLAGS
  to change your C linking flags
FC
  to choose your Fortran compiler
XXX_DIR
  CMake FindXXX will try this directory as a priority
XXX_ROOT
  CMake FindXXX will include this directory in the search

.. _modules: https://www.nersc.gov/users/software/user-environment/modules/

Tuning the configuration: ccmake
--------------------------------

When the configuration is successful, you can tune it using ccmake:

.. code: shell
  ccmake .

(notice the double c of ``ccmake``). This is an interactive tool, that lets you
choose the compilation parameters. Navigate with the arrows to the parameter you
want to change and hit enter to edit. Remember that any changes will be lost
when you invoke again a ``configure`` script.

Notable parameters are::

  BLA_VENDOR                      ALL (Typically you want either Intel10_64lp_seq, IBMESSL, or OpenBLAS)

Available in submodule PaRSEC builds only::

  PARSEC_DEBUG                    OFF (and all other PARSEC_DEBUG options)
  PARSEC_DIST_COLLECTIVES         ON
  PARSEC_GPU_WITH_CUDA            ON
  PARSEC_PROF_*                   OFF (all PARSEC_PROF_ flags off)

Using the *expert* mode (key 't' to toggle to expert mode), you can change other
useful options, like::

  CMAKE_C_FLAGS_RELEASE
  CMAKE_EXE_LINKER_FLAGS_RELEASE
  CMAKE_Fortran_FLAGS_RELEASE
  CMAKE_VERBOSE_MAKEFILE

And others to change the path to some compilers, for example. The
``CMAKE_VERBOSE_MAKEFILE`` option, when turned ``ON``, will display the command run when
compiling, which can help debugging configuration mistakes.  When you have set
all the options you want in ccmake, type 'c' to configure again, and 'g' to
generate the files. If you entered wrong values in some fields, ccmake will
complain at 'c' time.

BLAS and LAPACKE
================

Choosing a BLAS
---------------

DPLASMA needs to have access to a BLAS implementation and LAPACKE_ (+TMG)
interface. It is recommended that you use a vendor supplied BLAS (e.g., 
Intel MKL_, IBM ESSL_, OpenBLAS_, etc.) rather than a generic option. 
Using the reference BLAS_ (or, to a lesser extent, ATLAS_) often result in 
poor performance.

In order to control which BLAS will be selected, you can either

1. Pass the --with-blas=xxx to the configure script (see above)
2. Set the BLA_VENDOR CMake variable (-DBLA_VENDOR=xxx)

Typical values for these options are ``Intel10_64lp_seq`` (Intel MKL), ``IBMESSL``,
``OpenBLAS``, etc. You can refer to the CMake FindBLAS_ documentation to discover
more options.

.. _MKL: https://software.intel.com/en-us/mkl
.. _ESSL: https://www.ibm.com/support/knowledgecenter/en/SSFHY8/essl_welcome.html
.. _OpenBLAS: https://www.openblas.net
.. _ATLAS: http://math-atlas.sourceforge.net
.. _BLAS: https://github.com/Reference-LAPACK/lapack
.. _FindBLAS: https://cmake.org/cmake/help/latest/module/FindBLAS.html

LAPACKE
-------

LAPACKE lets C programs call Fortran LAPACK functions. Fortunately, many 
modern BLAS vendors (e.g., MKL, OpenBLAS) provide a full LAPACKE stack (including
CBLAS). In this case, just providing a BLAS is sufficient.

However, some vendors provide only a subset of LAPACK/LAPACKE (e.g., ESSL). In this
case, it is still recommended that you use the vendor BLAS, but you will need to
complement the missing features with the reference LAPACK/LAPACKE library.

.. code:: shell

  LAPACKE_ROOT=$LAPACK_BUILD_DIR configure --with-blas=IBMESSL 

OpenMP vs Serial BLAS
---------------------

In general, DPLASMA operates faster when using a serial BLAS, letting PaRSEC
manage parallelism. This setup can be achieved by linking with a serial version
of the BLAS library (``Intel10_64lp_seq`` rather than ``Intel10_64lp``), or 
alternatively, by disabling the OpenMP based BLAS-internal parallelism found in
many BLAS by setting the environment variable ``export OMP_NUM_THREADS=1`` at 
runtime.

Still, some architectures may benefit greatly from using an OpenMP BLAS, notably,
Intel KNC Phi accelerators on which OpenMP parallelism should be set to the number
of hardware threads per core. If you have an unusual architecture, experiment for
yourself!

Building DPLASMA
================

If the configuration was good, compilation should be as simple and
fancy as ``make``. To debug issues, use ``make VERBOSE=1`` or turn the
``CMAKE_VERBOSE_MAKEFILE`` option to ``ON`` using ``ccmake``. Check
your compilation lines, and adapt your configuration options accordingly.

Spack
-----

Some DOE sites are exploring the use of Spack_ to install software. You
can integrate PaRSEC in a Spack environment by using the provided
configurations in ``contrib/spack``. See the Readme there for more details.

Running with DPLASMA
====================

The dplasma library is compiled into ``dplasma/lib``. All testing programs are
compiled in ``dplasma/tests``. Examples are:

``dplasma/testing/testing_?getrf``
    LU Factorization (simple or double precision)
``dplasma/testing/testing_?geqrf``
    QR Factorization (simple or double precision)
``dplasma/testing/testing_?potrf``
    Cholesky Factorization (simple or double precision)

All the binaries should accept as input:

    -c <n>                  the number of threads used for kernel execution on each node.
                            This should be set to the number of cores. Remember that one
                            additional thread will be spawned to handle the communications
                            in the MPI version.
    -N SIZE                 a mandatory argument to define the size of the matrix
    -g <number of GPUs>     number of GPUs to use, if the operation is GPU-enabled
    -t <blocksize>          columns in a tile
    -T <blocksize>          rows in a tile, (WARNING: most algorithm included in DPLASMA
                            requires square tiles)
    -p <number of rows>     to require a 2-D block cyclic distribution of p rows
    -q <number of columns>  to require a 2D block cyclic distribution of q columns

A typical dplasma run using MPI looks like

.. code:: bash

  mpiexec -np 8 ./testing_spotrf -c 8 -g 0 -p 4 -q 2 -t 120 -T 120 -N 1000

This invocation run a Cholesky factorization on 8 nodes, 8 computing threads per node, nodes being
arranged in a ``4x2`` grid, with a distributed generation of the matrix of size ``1000x1000`` floats, with
tiles of size ``120x120``. Each test can dump the list of options with ``-h``. Some tests have specific options 
(like ``-I`` to tune the inner block size in QR and LU, and ``-M`` in LU or QR to have non-square matrices).

In addition to the parameters usually accepted by DPLASMA (see ``mpirun -np 1 ./testing_dpotrf --help`` for a full
list), the PaRSEC runtime engine can be tuned through its MCA. MCA parameters can be passed to the runtime engine
after the DPLASMA arguments, by separating the DPLASMA arguments from the PaRSEC arguments with -- (e.g. 
``mpirun -np 8 ./testing_dpotrf -c 8 -N 1000 -- --mca mca_sched ap`` would tell DPLASMA to use 8 cores, and PaRSEC 
to use the AP (Absolute Priority) scheduling heuristic). A complete list of MCA parameters can be found by passing 
``--help`` to the PaRSEC runtime engine (e.g. ``mpirun -np 1 ./testing_dpotrf -c 1 -N 100 -- --help``).

TL;DR
=====

.. code:: bash

  mkdir builddir && cd builddir
  ${srcdir}/configure --with-hwloc --with-mpi --with-blas=Intel10_64lp_seq --disable-debug --prefix=$PWD/install
  make install
  mpiexec -n 8 tests/testing_dpotrf -N 1000 -x -v

______

--
Happy hacking,
  The DPLASMA team.

