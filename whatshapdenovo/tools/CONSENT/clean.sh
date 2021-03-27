#!/bin/bash

cd BMEAN
cd Complete-Striped-Smith-Waterman-Library/src
make clean
cd ../../
cd spoa/build/
make clean
rm CMakeCache.txt
cd ../../
make clean
cd ../minimap2/
make clean
cd ../
make clean

