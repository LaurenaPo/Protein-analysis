from parse_stride import *
from parse_dssp import *
from parse_dssp2 import *
from missing_residues import *
from parse_sst import *
from math import log10
n = lambda x : int(log10(x))+1

def common_output(pdb_file, dssp_file, dssp2_file, stride_file,sst_file):
    common_file_parsed = dssp_file[0:4] + "_final_output.fasta"
    aa_pos, chain_ID_seq, aa_seq, dssp_struct_seq = parse_dssp(dssp_file)
    print(aa_pos)
    aa_pos_stride,stride_struct_seq = parse_stride(stride_file)
    dssp2_struct_seq= parse_dssp_pp2(dssp2_file)
    sst_seq = parse_promotif_sst(sst_file)
    missing_seq= missing_residues(pdb_file)
    missing_seq_1 = [int(i) for i in missing_seq]
    missing_seq_1=[i-1 for i in missing_seq_1]
    print(len(stride_struct_seq))
    print(len(aa_pos))
    with open(common_file_parsed, "w+") as file_out:
        file_out.write("  POS C A D 2 S P\n")
        x=0
        for i in range(len(stride_struct_seq)+len(missing_seq)):
            k=i-x
            if i in missing_seq_1:
                file_out.write(" {} {} {} {} {} {} {}\n".format(missing_seq[x], "$","$", "$",
                                                         "$", "$","$"))
                x+=1
            else:
                file_out.write(" {} {} {} {} {} {} {}\n".format(aa_pos[k], chain_ID_seq[k], aa_seq[k], dssp_struct_seq[k],
                                                         dssp2_struct_seq[k], stride_struct_seq[k], sst_seq[k]))



def main():
    print("1")
    pdb_file = "1oip.pdb"
    dssp_file = "1oip.dssp"
    dssp2_file = "1oip.dssp2"
    stride_file = "1oip.stride"
    sst_file = "1oip.sst"
    common_output(pdb_file, dssp_file,dssp2_file, stride_file, sst_file)


main()
exit(0)
