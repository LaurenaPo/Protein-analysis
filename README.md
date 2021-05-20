# Protein-analysis

This is a tool which gathers three tools to get an overview of the protein secondary structure assignment
The three tools are : DSSP, Stride and Promotif

## How to use
### Requirements
- DSSP, Stride and Promotif must been installed
- To launch the tool /bin/tcsh shell is needed
- Download the pdb file

To launch the script file. Use **./script.sh "pdb code"**  . This will launch the tool and create the output files.

## Output format
The main output file is named **"pdb code"_final_output** 
The file contains 10 columns

1. column 1 :
Position of the residue 
2. column 2 :
Protein chain identifier
3. column 3 :
One letter residue code
4. column 4 :
DSSP one letter secondary structure code
5. column 5 :
DSSP-PPII one letter secondary structure code
6. column 6 :
Stride one letter secondary structure code
7. column 7 :
Promotif one letter secondary structure code
8. column 8 :
Classical 3 states "H" "E" "C"
   * if the 3 tools give the same secondary structure code and match a classical 3 states code, the state will be in Capital letter
   * if only 2 tools give the same secondary structure code and match a classical 3 states code, the state will be in small letter
   * if only 1 tool match a classical 3 states code, the state will be :  ___*___
   * if no tools match a classical 3 states code, the state will be : ___-___
9. column 9 :
  A flag **#** can occur if between the 3 tools a difference is seen the difference can be :
   * difference of residue 
   * difference of secondary structure assignment 
   
* If a flag appears, see remarks at the end of the file

10.  column 10 :
Special structure from Promotif :
     * if the letter is **B** the residue form a beta-turn
     * if the letter is **G** the residue form a gamma-turn
     * if the letter is **L** the residue form a beta-bulge
     * if the letter is **M** the residue have multiple special structure

After the columns there is two confusion matrix. 
The first one is between DSSP and Stride
The second one is between DSSP and Promotif
