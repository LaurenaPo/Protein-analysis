from parse_stride import *
from parse_dssp import *
from parse_dssp2 import *
from missing_residues import *
from parse_sst import *
from parse_beta_turn import *
from math import log10
n = lambda x : int(log10(x))+1

def common_output(pdb_file, dssp_file, dssp2_file, stride_file,sst_file):
    common_file_parsed = dssp_file[0:4] + "_final_output.fasta"

    #recupere les donnee dssp
    aa_pos, chain_ID_seq, aa_seq, dssp_struct_seq = parse_dssp(dssp_file)

    #recupere les donnee stride
    aa_pos_stride, aa_seq_stride,stride_struct_seq = parse_stride(stride_file)

    #recupere les donnee dssp2
    dssp2_struct_seq= parse_dssp_pp2(dssp2_file)

    #recupere les donnee promotif_sst
    aa_pos_promotif, aa_seq_promotif,promotif_struct_seq= parse_promotif_sst(sst_file)

    #recupere les missing residues de la pdb
    missing_seq= missing_residues(pdb_file)
    missing_seq_1 = [int(i) for i in missing_seq]
    missing_seq_1=[i-1 for i in missing_seq_1]
    print(len(stride_struct_seq))
    print(len(aa_pos))

    #liste de position des flags
    flags=[]
    Remarque=""
    #compare dssp et stride
    compare_dssp_stride=compare_seq(aa_pos,aa_seq,aa_pos_stride, aa_seq_stride)
    print(compare_dssp_stride)
    if(len(compare_dssp_stride)>0):
        flags+=compare_dssp_stride

    #compare dssp et promotif
    compare_dssp_promotif=compare_seq(aa_pos,aa_seq,aa_pos_promotif, aa_seq_promotif)

    if(len(compare_dssp_promotif)>0):
        flags+=compare_dssp_promotif

    print(flags)
    flags=[int(i) for i in flags]
    print(flags)
    with open(common_file_parsed, "w+") as file_out:
        file_out.write("  POS C A D 2 S P\n")
        x=0
        j=0
        for i in range(len(stride_struct_seq)+len(missing_seq)):
            k=i-x
            if i in missing_seq_1:
                file_out.write(" {} {} {} {} {} {} {}\n".format(missing_seq[x], "$","$", "$",
                                                         "$", "$","$"))
                x+=1
            else:
                if k in flags:
                    print("1")
                    file_out.write(" {} {} {} {} {} {} {} {}\n".format(aa_pos[k], chain_ID_seq[k], aa_seq[k], dssp_struct_seq[k],
                                                         dssp2_struct_seq[k], stride_struct_seq[k], promotif_struct_seq[k],"#"))
                if k in compare_dssp_stride :
                    Remarque="#1 pos:{} chain:{} DSSP:{}  Stride:{}\n ".format(aa_pos[k],chain_ID_seq[k],aa_seq[k],aa_seq_stride[k])
                if k in compare_dssp_promotif:
                    Remarque="#2 pos:{} chain:{} DSSP:{}  Promotif:{}\n ".format(aa_pos[k],chain_ID_seq[k],aa_seq[k],aa_seq_promotif[k])
                file_out.write(" {} {} {} {} {} {} {}\n".format(aa_pos[k], chain_ID_seq[k], aa_seq[k], dssp_struct_seq[k],
                                                         dssp2_struct_seq[k], stride_struct_seq[k], promotif_struct_seq[k]))
        file_out.write(Remarque)



def main():
    pdb_file = "1oip.pdb"
    dssp_file = "1oip.dssp"
    dssp2_file = "1oip.dssp2"
    stride_file = "1oip.stride"
    sst_file = "1oip.sst"
    common_output(pdb_file, dssp_file,dssp2_file, stride_file, sst_file)


main()
exit(0)
