from parse_stride import *
from parse_dssp import *
from parse_dssp2 import *
from missing_residues import *
from parse_sst import *
from parse_beta_turn import *
from parse_gturn import *
from parse_bbulge import *
from classical_3_states import *
from math import log10
import sys

n = lambda x: int(log10(x)) + 1


def common_output(pdb_file, dssp_file, dssp2_file, stride_file, sst_file, bturn_file, gturn_file, bbulge_file):
    common_file_parsed = dssp_file[0:4] + "_final_output.txt"

    # recupere les donnee dssp
    aa_pos, chain_ID_seq, aa_seq, dssp_struct_seq = parse_dssp(dssp_file)

    # recupere les donnee stride
    aa_pos_stride, aa_seq_stride, stride_struct_seq = parse_stride(stride_file)

    # recupere les donnee dssp2
    dssp2_struct_seq = parse_dssp_pp2(dssp2_file)

    # recupere les donnee promotif_sst
    aa_pos_promotif, aa_seq_promotif, promotif_struct_seq = parse_promotif_sst(sst_file)

    # recupere les missing residues de la pdb
    missing_seq = missing_residues(pdb_file)
    missing_seq_1 = [int(i) for i in missing_seq]
    missing_seq_1 = [i - 1 for i in missing_seq_1]
    #print(len(stride_struct_seq))
    #print(len(aa_pos))

    # get classical 3 states
    classical_3_states_tab = classical_3_states(dssp_struct_seq, stride_struct_seq, promotif_struct_seq,
                                                len(dssp_struct_seq))
    #print("len classical 3 states = ", len(classical_3_states_tab))
    #print(classical_3_states_tab)

    # liste de position des flags
    flags = []
    other_struct = []
    Remarque1 = ""
    Remarque2 = ""
    Remarque3 = ""
    Remarque4 = ""
    Remarque5 = ""
    Remarque6 = ""
    Remarque7 = ""
    # compare dssp et stride acides amine
    compare_dssp_stride = compare_seq(aa_pos, aa_seq, aa_pos_stride, aa_seq_stride)

    if (len(compare_dssp_stride) > 0):
        flags += compare_dssp_stride

    # compare dssp et promotif acides amine
    compare_dssp_promotif = compare_seq(aa_pos, aa_seq, aa_pos_promotif, aa_seq_promotif)
    if (len(compare_dssp_promotif) > 0):
        flags += compare_dssp_promotif

    # compare dssp et stride structure secondaire
    compare_dssp_stride_struct = compare_seq(aa_pos, dssp_struct_seq, aa_pos_stride, stride_struct_seq)
    if (len(compare_dssp_stride_struct) > 0):
        flags += compare_dssp_stride_struct

    # compare dssp et stride structure secondaire
    compare_dssp_promotif_struct = compare_seq(aa_pos, dssp_struct_seq, aa_pos_promotif, promotif_struct_seq)
    if (len(compare_dssp_promotif_struct) > 0):
        flags += compare_dssp_promotif_struct

    # matrice de confusion dssp promotif
    matrice_dssp_promotif = matrice_seq_8etats(dssp_struct_seq, promotif_struct_seq)
    affichage_matrice_dssp_promotif = "              Matrice de confusion DSSP et Prmotif\n" + print_matrix_8etats(
        matrice_dssp_promotif)

    # matrice de confusion dssp promotif
    matrice_dssp_stride = matrice_seq_7etats(dssp_struct_seq, stride_struct_seq)
    affichage_matrice_dssp_stride = "              Matrice de confusion DSSP et Stride\n" + print_matrix_7etats(
        matrice_dssp_stride)

    # beta-turn
    bturn_sequence, aa_pos_bturn, bturn_type, chain_ID_bturn = parse_promotif_bturn(bturn_file)
    aa_pos_bturn = [int(i) for i in aa_pos_bturn]

    if (len(aa_pos_bturn) > 0):
        other_struct += aa_pos_bturn

    # gamma-turn
    gturn_sequence, aa_pos_gturn, gturn_type, chain_ID_gturn = parse_promotif_gturn(gturn_file)
    aa_pos_gturn = [int(i) for i in aa_pos_gturn]

    if (len(aa_pos_gturn) > 0):
        other_struct += aa_pos_gturn

    # beta-bulges
    x_seq_bbulge, x_pos_bbulge, first_seq_bbulge, first_pos_bbulge, second_seq_bbulge, second_pos_bbulge, bbluge_type, chain_ID_bbulge = parse_promotif_bbulge(
        bbulge_file)
    x_pos_bbulge = [int(i) for i in x_pos_bbulge]
    first_pos_bbulge = [int(i) for i in first_pos_bbulge]
    second_pos_bbulge = [int(i) for i in second_pos_bbulge]

    if (len(x_pos_bbulge) > 0):
        other_struct += x_pos_bbulge
        other_struct += first_pos_bbulge
        other_struct += second_pos_bbulge

    # cherche les positions qui ont deux assignations structurales speciales
    once = set()
    seenOnce = once.add
    twice = set(num for num in other_struct if num in once or seenOnce(num))
    twice = list(twice)

    flags = [int(i) for i in flags]
    #print(flags)
    with open(common_file_parsed, "w+") as file_out:
        file_out.write("  POS C A D 2 S P\n")
        x = 0
        for i in range(len(stride_struct_seq) + len(missing_seq)):
            k = i - x

            if i in missing_seq_1:
                file_out.write(" {} {} {} {} {} {} {}\n".format(missing_seq[x], "$", "$", "$",
                                                                "$", "$", "$"))
                x += 1
            else:
                if k in flags:

                    file_out.write(
                        " {} {} {} {} {} {} {} {} {}".format(aa_pos[k], chain_ID_seq[k], aa_seq[k], dssp_struct_seq[k],
                                                             dssp2_struct_seq[k], stride_struct_seq[k],
                                                             promotif_struct_seq[k], classical_3_states_tab[k], "#"))

                    if i + 1 in aa_pos_bturn:
                        if i + 1 in twice:
                            file_out.write(" M\n")
                        else:
                            file_out.write(" B\n")
                    elif i + 1 in aa_pos_gturn:
                        if i + 1 in twice:
                            file_out.write(" M\n")
                        else:
                            file_out.write(" G\n")
                    elif i + 1 in x_pos_bbulge:
                        if i + 1 in twice:
                            file_out.write(" M\n")
                        else:
                            file_out.write(" L\n")
                    elif i + 1 in first_pos_bbulge:
                        if i + 1 in twice:
                            file_out.write(" M\n")
                        else:
                            file_out.write(" L\n")
                    elif i + 1 in second_pos_bbulge:
                        if i + 1 in twice:
                            file_out.write(" M\n")
                        else:
                            file_out.write(" L\n")
                    else:
                        file_out.write("\n")

                    if k in compare_dssp_stride:
                        Remarque1 = Remarque1 + "Remark #1 pos:{} chain:{} DSSP:{}  Stride:{}\n".format(aa_pos[k],
                                                                                                        chain_ID_seq[k],
                                                                                                        aa_seq[k],
                                                                                                        aa_seq_stride[
                                                                                                            k])

                    if k in compare_dssp_promotif:
                        Remarque2 = Remarque2 + "Remark #2 pos:{} chain:{} DSSP:{}  Promotif:{}\n".format(aa_pos[k],
                                                                                                          chain_ID_seq[
                                                                                                              k],
                                                                                                          aa_seq[k],
                                                                                                          aa_seq_promotif[
                                                                                                              k])

                    if k in compare_dssp_stride_struct:
                        Remarque3 = Remarque3 + "Remark #3 pos:{} chain:{} DSSP:{}  Stride:{}\n".format(aa_pos[k],
                                                                                                        chain_ID_seq[k],
                                                                                                        dssp_struct_seq[
                                                                                                            k],
                                                                                                        stride_struct_seq[
                                                                                                            k])

                    if k in compare_dssp_promotif_struct:
                        Remarque4 = Remarque4 + "Remark #4 pos:{} chain:{} DSSP:{}  Promotif:{}\n".format(aa_pos[k],
                                                                                                          chain_ID_seq[
                                                                                                              k],
                                                                                                          dssp_struct_seq[
                                                                                                              k],
                                                                                                          promotif_struct_seq[
                                                                                                              k])

                else:
                    file_out.write(
                        " {} {} {} {} {} {} {} {}".format(aa_pos[k], chain_ID_seq[k], aa_seq[k], dssp_struct_seq[k],
                                                          dssp2_struct_seq[k], stride_struct_seq[k],
                                                          promotif_struct_seq[k], classical_3_states_tab[k]))
                    if i + 1 in aa_pos_bturn:
                        if i + 1 in twice:
                            file_out.write("   M\n")
                        else:
                            file_out.write("   B\n")
                    elif i + 1 in aa_pos_gturn:
                        if i + 1 in twice:
                            file_out.write("   M\n")
                        else:
                            file_out.write("   G\n")
                    elif i + 1 in x_pos_bbulge:
                        if i + 1 in twice:
                            file_out.write("   M\n")
                        else:
                            file_out.write("   L\n")
                    elif i + 1 in first_pos_bbulge:
                        if i + 1 in twice:
                            file_out.write("   M\n")
                        else:
                            file_out.write("   L\n")
                    elif i + 1 in second_pos_bbulge:
                        if i + 1 in twice:
                            file_out.write("   M\n")
                        else:
                            file_out.write("   L\n")
                    else:
                        file_out.write("\n")

        for i in range(len(aa_pos_bturn) + 1):
            if (i % 4 == 0 and i != 0):
                Remarque5 = Remarque5 + "Remark #5 beta-turn chain : {} pos : {}-{}-{}-{} sequence : {} type: {}\n".format(
                    chain_ID_bturn[(i // 4) - 1], aa_pos_bturn[i - 4], aa_pos_bturn[i - 3], aa_pos_bturn[i - 2],
                    aa_pos_bturn[i - 1],
                    bturn_sequence[(i // 4) - 1], bturn_type[(i // 4) - 1])

        for i in range(len(aa_pos_gturn) + 1):
            if (i % 3 == 0 and i != 0):
                Remarque6 = Remarque6 + "Remark #6 gamma-turn chain : {} pos : {}-{}-{} sequence : {} type: {}\n".format(
                    chain_ID_gturn[(i // 3) - 1], aa_pos_gturn[i - 3], aa_pos_gturn[i - 2], aa_pos_gturn[i - 1],
                    gturn_sequence[(i // 3) - 1], gturn_type[(i // 3) - 1])

        for i in range(len(x_pos_bbulge)):
            Remarque7 = Remarque7 + "Remark #7 beta-bulges chain : {} pos-x: {} pos1-2 : {}-{} sequence-x : {} sequence1-2 : {}-{} type: {}\n".format(
                chain_ID_gturn[i], x_pos_bbulge[i], first_pos_bbulge[i], second_pos_bbulge[i], x_seq_bbulge[i],
                first_seq_bbulge[i], second_seq_bbulge[i], bbluge_type[i])

        file_out.write("\n" + affichage_matrice_dssp_stride + "\n")
        file_out.write(affichage_matrice_dssp_promotif + "\n")
        file_out.write(Remarque1 + "\n")
        file_out.write(Remarque2 + "\n")
        file_out.write(Remarque3 + "\n")
        file_out.write(Remarque4 + "\n")
        file_out.write(Remarque5 + "\n")
        file_out.write(Remarque6 + "\n")
        file_out.write(Remarque7 + "\n")


def main():
    args = sys.argv[0:]
    if len(args[1])==4:
        pdb_file = args[1]+".pdb"
        dssp_file = pdb_file[0:4] + ".dssp"
        dssp2_file = pdb_file[0:4] + ".dssp2"
        stride_file = pdb_file[0:4] + ".stride"
        sst_file = pdb_file[0:4] + ".sst"
        bturn_file = pdb_file[0:4] + ".bturns"
        gturn_file = pdb_file[0:4] + ".gturns"
        bbulge_file = pdb_file[0:4] + ".blg"
        common_output(pdb_file, dssp_file, dssp2_file, stride_file, sst_file, bturn_file, gturn_file, bbulge_file) 
    else:
        print("fichier {}.pdb incorrect".format(args[1]))


main()
exit(0)
