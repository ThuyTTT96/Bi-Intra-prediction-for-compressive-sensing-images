#!/bin/bash


CURRENT_PATH="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
BUILD_PATH="../../libs/KLab/main/macos/"


cd $BUILD_PATH
. UnixMakefile_x86Debug_Make.sh

cd $CURRENT_PATH
cd $BUILD_PATH
cd UnixMakefile_x86Debug
. build.sh

cd $CURRENT_PATH
cd $BUILD_PATH
. UnixMakefile_x86Release_Make.sh

cd $CURRENT_PATH
cd $BUILD_PATH
cd UnixMakefile_x86Release
. build.sh
