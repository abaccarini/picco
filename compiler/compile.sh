#!/bin/bash

mkdir -p bin # making bin directory if does not exist
cd src/utility
make
mv picco-utility ../../bin/
