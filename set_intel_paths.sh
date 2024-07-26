
################################################################################
# The following commands create directories for Intel's compiler suite

# Intel compiler installation directory
export intel_root=$HOME/intel/oneapi

# Intel compiler directory at $HOME
export intel_dir=$HOME/intel

# Make directory
rm -f $intel_dir/include $intel_dir/lib
rm -r $intel_dir/bin $intel_dir/mkl $intel_dir/mpi
mkdir -p $intel_dir
mkdir -p $intel_dir/bin
mkdir -p $intel_dir/mkl
mkdir -p $intel_dir/mpi

# Make links
cd $intel_dir
ln -sfn $intel_root/compiler/latest/lib/ $intel_dir/lib

ln -sfn $intel_root/compiler/latest/bin/ifort $intel_dir/bin/ifort
ln -sfn $intel_root/compiler/latest/bin/icc $intel_dir/bin/icc
ln -sfn $intel_root/compiler/latest/bin/icpc $intel_dir/bin/icpc
ln -sfn $intel_root/compiler/latest/bin/ifx $intel_dir/bin/ifx
ln -sfn $intel_root/compiler/latest/bin/icx $intel_dir/bin/icx
ln -sfn $intel_root/compiler/latest/bin/icpx $intel_dir/bin/icpx

ln -sfn $intel_root/mkl/latest/include/ $intel_dir/mkl/include
ln -sfn $intel_root/mkl/latest/lib/intel64/ $intel_dir/mkl/lib

ln -sfn $intel_root/mpi/latest/bin/ $intel_dir/mpi/bin
ln -sfn $intel_root/mpi/latest/include/mpi $intel_dir/mpi/include
ln -sfn $intel_root/mpi/latest/lib/ $intel_dir/mpi/lib
ln -sfn $intel_root/mpi/latest/env/ $intel_dir/mpi/env

################################################################################
# Due to changes in Intel's compiler suite, the following paths are no longr needed
#mkdir -p $intel_dir/include
#ln -sfn $intel_root/compiler/latest/include/ $intel_dir/include

exit

################################################################################
# The following system paths should be set. 
# Depending on your shell environment, you may need to adjust how to set these paths. 
# The following commands are for bash.

export PATH=$PATH:$HOME/intel/bin:$HOME/intel/mpi/bin
export PATH=$PATH:$HOME/intel/mkl/include:$HOME/intel/mpi/include
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/intel/lib:$HOME/intel/mkl/lib:$HOME/intel/mpi/lib:$HOME/intel/mpi/lib/release

# $HOME/intel/include
