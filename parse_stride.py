from ressources import *

def parse_stride(stride_file):
    #Recovery of protein and structural sequences in .stride file
    aa_pos=[]
    aa_seq = ""
    struct_seq = ""
    chain_ID_seq = ""
    with open(stride_file, 'r') as file_in:
        for line in file_in:
            if line[0:3] == "ASG":
                aa_pos.append(line[11:15])
                aa_seq+=AA[line[5:8]]
                struct_seq+=line[24:25]
                chain_ID_seq += line[9:10]

    #Writing sequences in a fasta file
    stride_file_parsed = stride_file[0:4]+"_output_stride.fasta"
    with open(stride_file_parsed, "w+") as file_out:
        file_out.write(">{} PARSED\n".format(stride_file))
        #Column1: position
        #Column2: amino acid sequence in one letter code
        #Column3: one-letter chain ID, if any
        #Column4: secondary structure
        file_out.write("  POS C A S\n")
        for i in range(len(aa_seq)):
            file_out.write(" {} {} {} {}\n".format(aa_pos[i], chain_ID_seq[i], aa_seq[i], struct_seq[i]))
    return aa_pos, aa_seq, struct_seq

"""def main():
    stride_file = "1oip.stride"
    aa_pos,aa_seq, struc_seq = parse_stride(stride_file)
    print(aa_pos)
    print(aa_seq)
    print(struc_seq)

main()"""
