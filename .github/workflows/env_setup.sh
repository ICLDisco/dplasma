# This file should be "sourced" into your environment

# Show the executed command, but don't affect spawned shells
trap 'echo "# $BASH_COMMAND"' DEBUG # Show commands

echo "Loading environment"
if [[ -z "$SPACK_SETUP" || ! -e "$SPACK_SETUP" ]]; then
   echo Error! Environment variable \$SPACK_SETUP must point
   echo to a valid setup-env.sh Spack setup script.
   exit 1
fi
source $SPACK_SETUP
spack env activate dplasma

HIP=OFF
CUDA=OFF
if [ "$DEVICE" = "gpu_nvidia" ]; then
   spack load cuda
   CUDA=ON
elif [ "$DEVICE" = "gpu_amd" ];then
   HIP=ON
fi


# Disable RECURSIVE in CI tests until a real solution to https://github.com/ICLDisco/parsec/issues/548 is implemented
! read -d '' BUILD_CONFIG << EOF
        -G Ninja
        -DCMAKE_BUILD_TYPE=$BUILD_TYPE
        -DBUILD_SHARED_LIBS=$SHARED_TYPE
        -DMPIEXEC_PREFLAGS='--bind-to;none;--oversubscribe'
        -DCMAKE_INSTALL_PREFIX=$INSTALL_DIRECTORY
        -DDPLASMA_PRECISIONS=d
        -DPARSEC_HAVE_DEV_RECURSIVE_SUPPORT=OFF
        -DDPLASMA_GPU_WITH_CUDA=$CUDA
        -DDPLASMA_GPU_WITH_HIP=$HIP
EOF
export BUILD_CONFIG
