#!/bin/bash

# Record commands
set -o xtrace

# Make POA visualizations of each difference
echo Making POA graph visualizations of difference
mkdir visualizations
for i in $(cat ${1} | grep 'Mismatch' | cut -d' ' -f4); do
  echo hello ${i}
  python3 ../../scripts/phasedPoaToDot.py ${2} --start=${i} --length=10 --output=visualizations/diff_${i}
done
