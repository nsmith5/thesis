#!/bin/bash
rm workspace/input workspace/benchmark workspace/init.yaml
cd lib
make
cd ../examples
make
cd ..
cp -v examples/input examples/benchmark init.yaml tools/h5topng.py workspace
