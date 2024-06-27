
################################################################################
# The following commands create directories for Intel's compiler suite

# Intel compiler installation directory
export intel_root=$HOME/intel/oneapi

# Intel compiler directory at $HOME, for convenience
export intel_dir=$HOME/intel

# Make directory
mkdir -p $intel_dir
mkdir -p $intel_dir/bin
mkdir -p $intel_dir/mkl
mkdir -p $intel_dir/mpi
mkdir -p $intel_dir/include

# Make links
cd $intel_dir
ln -sfn $intel_root/compiler/latest/linux/compiler/lib/intel64/ ./lib

ln -sfn $intel_root/compiler/latest/linux/compiler/include/ ./include

ln -sfn $intel_root/compiler/latest/linux/bin/intel64/ifort ./bin/ifort
ln -sfn $intel_root/compiler/latest/linux/bin/intel64/icc ./bin/icc
ln -sfn $intel_root/compiler/latest/linux/bin/intel64/icpc ./bin/icpc
ln -sfn $intel_root/compiler/latest/linux/bin/ifx ./bin/ifx
ln -sfn $intel_root/compiler/latest/linux/bin/icx ./bin/icx
ln -sfn $intel_root/compiler/latest/linux/bin/icpx ./bin/icpx

ln -sfn $intel_root/mkl/latest/include/ ./mkl/include
ln -sfn $intel_root/mkl/latest/lib/intel64/ ./mkl/lib

ln -sfn $intel_root/mpi/latest/bin/ ./mpi/bin
ln -sfn $intel_root/mpi/latest/include/ ./mpi/include
ln -sfn $intel_root/mpi/latest/lib/ ./mpi/lib
ln -sfn $intel_root/mpi/latest/env/ ./mpi/env

################################################################################
# The following system paths should be set. 
# Depending on your shell environment, you may need to adjust how to set these paths. 
# The following commands are for bash.

exit

export PATH=$PATH:$HOME/intel/bin:$HOME/intel/mpi/bin
export PATH=$PATH:$HOME/intel/include:$HOME/intel/mkl/include:$HOME/intel/mpi/include
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/intel/lib:$HOME/intel/mkl/lib:$HOME/intel/mpi/lib:$HOME/intel/mpi/lib/release
