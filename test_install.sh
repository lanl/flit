
# Set Intel's coompiler paths and patch paths to system shell configuration 
./set_intel_paths.sh --intel-base $HOME/tool/intel/oneapi --intel-root $HOME/intel --patch-shell 

# **************** Shell restart is required here ************************

# Install FLIT with HDF5
ruby install.rb --intel-root $HOME/intel --with-hdf5 --hdf5-root $HOME/tool/hdf5 

# # Rebuild FLIT with HDF5 
# ruby install.rb --intel-root $HOME/intel --with-hdf5 --hdf5-root $HOME/tool/hdf5 --clean

# # Install FLIT with more specifications
# ruby install.rb --intel-root $HOME/intel --with-hdf5 --hdf5-root $HOME/tool/hdf5 --hdf5-version 2.1.0 --build-root /tmp/h5

# # Install FLIT without HDF5
# ruby install.rb --intel-root $HOME/intel --without-hdf5
