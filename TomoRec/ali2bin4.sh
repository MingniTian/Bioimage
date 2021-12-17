#!/bin/bash
#
for f in ./*.ali; do
prfx=`basename ${f} .ali`
  echo; echo ${prfx}
  binvol -binning 4 -zbinning 1 ${prfx}.ali ${prfx}.ali.bin4
done
