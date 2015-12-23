#!/bin/bash

# Checks that user provided the required input.
if [ -z "$1" -o -z "$2" ]; then
  echo 'Usage: ./make_system [# blocks/chain] [# chains]'
  exit 1
fi

python bead_spring_system.py  \
   --num_chains $2 \
   --block SSSSSSCAUBBUACSSSSSSS \
   --num-block $1 \
   --filename $1blk-$2chain.lammps \
   --bond_length "AC=3.628750,AU=3.804510,BB=4.894340,BU=3.811221,CS=4.904590,SS=4.959297" \
   --angle "ACS=135.011397,AUB=117.388454,BBU=149.864226,CAU=162.244489,CSS=139.965775,SSS=143.976966" \
   --density 1.0713 \
   --bead_masses "A=76.09364,U=58.03989,B=83.112325,C=72.0593,S=72.10776" \
   --seed 1234



