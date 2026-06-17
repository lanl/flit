
# System shell configuration file
shell_config=~/.bashrc

# Intel base directory
intel_base=$HOME/tool/intel/oneapi

# Set Intel's compiler paths, patch the selected shell configuration
./set_intel_paths.sh --intel-base $intel_base --intel-root $HOME/intel --patch-shell --shell-config $shell_config

# Update shell configuration or just restart the shell
source $shell_config

# Install FLIT with HDF5
ruby install.rb --intel-root $HOME/intel --with-hdf5 --hdf5-root $HOME/tool/hdf5 --shell-config $shell_config

# # Rebuild FLIT with HDF5 
# ruby install.rb --intel-root $HOME/intel --with-hdf5 --hdf5-root $HOME/tool/hdf5 --shell-config $shell_config --clean

# # Install FLIT with more specifications
# ruby install.rb --intel-root $HOME/intel --with-hdf5 --hdf5-root $HOME/tool/hdf5 --hdf5-version 2.1.0 --build-root /tmp/h5 --shell-config $shell_config

# # Install FLIT without HDF5
# ruby install.rb --intel-root $HOME/intel --without-hdf5

# Update shell configuration or just restart the shell
source $shell_config
