#!/bin/bash

rm -rf bin
mkdir -p bin
cd BMEAN;
./install.sh;
cd ..;
make;
