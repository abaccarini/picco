#!/bin/bash

mode="-DDEPLOYMENT=ON"
technique="-DSHAMIR=ON"

mkdir -p build
cd build

make clean
cmake -DCMAKE_BUILD_TYPE=Release $mode $technique ..
make -j24
cd ..
