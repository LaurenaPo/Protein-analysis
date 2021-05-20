#!/bin/tcsh
set file=$1
echo ${file}
dssp -i ${file}.pdb -o ${file}.dssp
perl dssppII_new.pl ${file}.pdb > ${file}.dssp2
stride ${file}.pdb > ${file}.stride
promotif ${file}.pdb
python3 ../final.py ${file}