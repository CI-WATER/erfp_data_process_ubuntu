#!/bin/sh
NAME="alan"
NOW=$(date +"%Y%m%d%H%M%S")
#RAPID environmental variables
export TACC_NETCDF_LIB='/home/$NAME/installz/netcdf-3.6.3-install/lib'
export TACC_NETCDF_INC='/home/$NAME/installz/netcdf-3.6.3-install/include'
export PETSC_DIR='/home/$NAME/installz/petsc-3.3-p7'
export PETSC_ARCH='linux-gcc-cxx'
#export PETSC_ARCH='linux-gcc-cxx-debug'
export TAO_DIR='/home/$NAME/installz/tao-2.1-p2'
export PATH=$PATH:/$PETSC_DIR/$PETSC_ARCH/bin
export PATH=$PATH:/home/$NAME/installz/netcdf-3.6.3-install/bin
#start tethys virtualenv
. /usr/lib/tethys/bin/activate
#run script from python
/usr/lib/tethys/bin/python /home/$NAME/work/scripts/erfp_data_process_ubuntu/rapid_process_async_ubuntu.py 1> /home/$NAME/work/logs/ecmwf_rapid_process_$NOW.log 2>&1
