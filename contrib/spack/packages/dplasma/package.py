##############################################################################
# Copyright (c) 2019       The University of Tennessee and the University
#                          of Tennessee Research Foundation.  All rights
#                          reserved.
#
# $COPYRIGHT$
#
##############################################################################
from spack import *

class Dplasma(CMakePackage):
    """DPLASMA: The distributed, accelerator enabled, dense linear algebra library"""

    homepage = "https://bitbucket.org/icldistcomp/dplasma"
    url      = "https://bitbucket.org/icldistcomp/dplasma/get/v2.0.0.tar.bz2"

    version('devel', git='https://bitbucket.org/icldistcomp/dplasma/git', branch='master')
    version('2.0.0', '')

    # We always need MPI for now.
    #variant('mpi', default=True, description='Use MPI for dependency transport between nodes')
    variant('cuda', default=True, description='Use CUDA for GPU acceleration')
    variant('profile', default=False, description='Generate profiling data')
    variant('debug', default=False, description='Debug version (incurs performance overhead!)')
    #variant('xcompile', default=False, description='Cross compile')

    depends_on('parsec')
    depends_on('cuda', when='+cuda')

    def configure_args(self):
        spec = self.spec
        return [
            '-DCMAKE_BUILD_TYPE=%s' % ('Debug' if '+debug' in spec else 'RelWithDebInfo')
#            '-DCUDA_TOOLKIT_ROOT_DIR=%s' %
#            '-DPARSEC_DIST_WITH_MPI=%s' % ('YES' if '-mpi' in spec else 'NO'),
#            '-DMPI_C_COMPILER=%s'
#            '-DMPI_CXX_COMPILER=%s'
#            '-DMPI_Fortran_COMPILER=%s'
        ]

