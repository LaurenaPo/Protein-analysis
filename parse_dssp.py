from ressources import *

def parse_dssp(dssp_file):
    #Recovery of protein and structural sequences in .dssp file

    #Initialization:
    aa_pos = []
    chain_ID_seq = ""
    aa_seq = ""
    struct_seq = ""
    remark_dict = {}
    phi_list = []
    psi_list = []

    #Reading of the file:
    with open(dssp_file, 'r') as file_in:
        #We are looking for the field # to start parsing
        line = file_in.readline().split()
        while line[0] != "#":
            line = file_in.readline().split()
        line = file_in.readline()

        #Parsing:
        while line != "":
            #original PDB resname:

            aa_code = line[13]
            if aa_code != "!":
                aa_pos.append(line[6:10])
                chain_ID_seq += line[11]

                if aa_code in aa_code_list:
                    aa_seq += aa_code
                elif aa_code == " ":
                    aa_seq += " "
                elif aa_code == "!":
                    aa_seq += " "
                else:
                    #There is a problem
                    aa_seq += "$"

                struct_code = line[16]
                if struct_code in struct_code_list:
                    struct_seq += struct_code
                elif struct_code == " ":
                    struct_seq += "C"
                else:
                    #There is a problem
                    struct_seq += "$"

                phi_list.append(line[103:109].strip())
                psi_list.append(line[109:115].strip())

            line = file_in.readline()

    #print("phi: ", phi_list)
    #print("psi: ", psi_list)
    #Writing sequences in a fasta file
    dssp_file_parsed = dssp_file[0:4] + "_output_dssp.txt"
    with open(dssp_file_parsed, "w+") as file_out:
        file_out.write(">{} PARSED\n".format(dssp_file))
        #Column1: position
        #Column2: amino acid sequence in one letter code
        #Column3: one-letter chain ID, if any
        #Column4: secondary structure
        file_out.write("  POS C A S\n")
        for i in range(len(aa_seq)):
            file_out.write(" {} {} {} {}\n".format(aa_pos[i], chain_ID_seq[i], aa_seq[i], struct_seq[i]))

    return aa_pos, chain_ID_seq, aa_seq, struct_seq

"""def main():
    dssp_file = "1oip.dssp"
    aa_pos, chain_ID_seq, aa_seq, struc_seq = parse_dssp(dssp_file)
    print(aa_pos)
    print(chain_ID_seq)
    print(aa_seq)
    print(struc_seq)

main()
exit(0)"""