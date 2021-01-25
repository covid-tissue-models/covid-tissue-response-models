#!/bin/bash

# Default argument values
PASSED_STATUS=0
GENERATE_STATUS=0

# Process arguments
for arg in "$@"
do
    case $arg in
        -r|--run_script)
        RUNSCRIPT="$2"
        echo $RUNSCRIPT
        shift # Remove name from processing
        shift # Remove value from processing
        ;;
        -i|--input)
        SCRIPTINPUT="$2"
        echo $SCRIPTINPUT
        shift # Remove name from processing
        shift # Remove value from processing
        ;;
        -s|--status)
        GENERATE_STATUS=0
        PASSED_STATUS=1
        STATUSFILE="$2"
        echo $STATUSFILE
        shift # Remove name from processing
        shift # Remove value from processing
        ;;
        -g|--generate_status)
        GENERATE_STATUS=1
        shift # Remove argument from processing
        ;;
    esac
done

# Set PREFIX_CC3D to the root directory of your CC3D installation!

#export PREFIX_CC3D=/geode2/home/u100/tjsego/Carbonate/cc3d/423_auto
export PREFIX_CC3D=/N/u/jferrari/Carbonate/cc3d_compile/vMaster/install

# If Python is not installed in the CC3D root directory, then replace this with the path to the Python root directory
# export PYTHON_INSTALL_PATH=${PREFIX_CC3D}/Python37/bin
#export PYTHON_INSTALL_PATH=/geode2/home/u100/tjsego/Carbonate/.conda/envs/cc3d_2020/bin
export PYTHON_INSTALL_PATH=/N/u/jferrari/Carbonate/.conda/envs/cc3d_master/bin


# Environment setup

export LC_NUMERIC="C.UTF-8"

export PATH=$PYTHON_INSTALL_PATH:$PATH
export LD_LIBRARY_PATH=${PREFIX_CC3D}/lib/:$LD_LIBRARY_PATH

export COMPUCELL3D_PLUGIN_PATH=${PREFIX_CC3D}/lib/site-packages/cc3d/cpp/CompuCell3DPlugins
export COMPUCELL3D_STEPPABLE_PATH=${PREFIX_CC3D}/lib/site-packages/cc3d/cpp/CompuCell3DSteppables

export LD_LIBRARY_PATH=${PREFIX_CC3D}/lib/site-packages/cc3d/cpp/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=${COMPUCELL3D_PLUGIN_PATH}:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=${COMPUCELL3D_STEPPABLE_PATH}:$LD_LIBRARY_PATH

export COMPUCELL3D_MAJOR_VERSION=4
export COMPUCELL3D_MINOR_VERSION=1
export COMPUCELL3D_BUILD_VERSION=0

export PYTHONPATH=${PREFIX_CC3D}/lib/site-packages

# Run it!
export exit_code=0
if [ "$GENERATE_STATUS" -eq 1 ] ; then
  # Run and generate status file
  ${PYTHON_INSTALL_PATH}/python $RUNSCRIPT -i $SCRIPTINPUT -g
elif [ "$PASSED_STATUS" -eq 1 ] ; then
  # Run and use passed status file
  ${PYTHON_INSTALL_PATH}/python $RUNSCRIPT -i $SCRIPTINPUT -s $STATUSFILE
else
  # Just run
  ${PYTHON_INSTALL_PATH}/python $RUNSCRIPT -i $SCRIPTINPUT
fi
exit_code=$?
