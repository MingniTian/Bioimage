#!/bin/bash
for f in ./*.ali.bin4; do

  prfx=`basename ${f} .ali.bin4`
  #echo; echo ${prfx}
  #tomo3d is used to generate 3D maps quickly
   mtffilter -l '0.20 0.05' ${prfx}.ali.bin4 ${prfx}.ali.low4
   tomo3d -a ${prfx}.tlt -i ${prfx}.ali.low4 -z 500 -S -o ${prfx}_tomo3d.low4 -m 0.20 
   clip rotx ${prfx}_tomo3d.low4 ${prfx}_tomo3d_xyz.low4
   rm -f ${prfx}_tomo3d.low4 ${prfx}.ali.low4

done
