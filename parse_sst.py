
def parse_promotif_sst(sst_file):
	aa_sequence = []
	aa_pos = []
	chain_ID = []
	struct_seq = []
	with open(sst_file, 'r') as file_in:
		word = "Sec"
		#Skip th wo first lines
		for i in range(7):
			next(file_in)
		for line in file_in:
			if word not in line:
				aa_pos.append(line[8:11])
				aa_sequence.append(line[19:20])
				chain_ID.append(line[6:7])
				if line[21:22] ==" ":
					struct_seq.append("C")
				else:
					struct_seq.append(line[21:22])
				print(" {} Chain->{} pos ->{}  {}  {} {}".format(line[1:4],line[6:7],line[7:11],line[12:19],line[19:20],line[21:22]))
			else:
				break
    #Writing sequences in a fasta file
	sst_parsed_file = sst_file[0:4]+"Output_sst.fasta"
	with open(sst_parsed_file,"w+") as file_out:
		file_out.write(">{} PARSED\n\n".format(sst_file))
		file_out.write("  Position   Chain   Sequence       Types\n")
		#column1: position
		#Column2: one-letter chain ID, if any
		#Column3: amino acid sequence in one letter code
		for i in range(len(aa_sequence)):
			file_out.write("   {}         {}       {}         {}\n".format(aa_pos[i], chain_ID[i], aa_sequence[i],struct_seq[i]))
	return struct_seq


"""def main():
	bturn_file = "1oip.sst"
	sequence = parse_promotif_sst(bturn_file)"""

		

#main()
#exit(0)
