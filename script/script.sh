#!/bin/tcsh
dssp -i 1oip.pdb -o 1oip.dssp;
perl dssppII_new.pl 1oip.pdb > 1oip.dssp2;
stride 1oip.pdb > 1oip.stride;
promotif 1oip.pdb;
python3 ../final.py;