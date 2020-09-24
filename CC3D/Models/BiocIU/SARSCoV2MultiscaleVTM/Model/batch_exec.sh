#!/bin/bash

# Set PREFIX_CC3D to the root directory of your CC3D installation!
export PREFIX_CC3D=/geode2/home/u100/tjsego/Carbonate/cc3d/423_auto

# If Python is not installed in the CC3D root directory, then replace this with the path to the Python root directory
# export PYTHON_INSTALL_PATH=${PREFIX_CC3D}/Python37/bin
export PYTHON_INSTALL_PATH=/geode2/home/u100/tjsego/Carbonate/.conda/envs/cc3d_2020/bin

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
BATCHSCRIPT=$(cd `dirname $0` && pwd)/batch_exec.py
echo $BATCHSCRIPT
${PYTHON_INSTALL_PATH}/python $BATCHSCRIPT
exit_code=$?
