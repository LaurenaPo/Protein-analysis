def parse_promotif_bturn(bturn_file):
	bturn_sequence = []
	aa_pos = []
	chain_ID = []
	bturn_type = []
	with open(bturn_file, 'r') as file_in:
		for i in range(2):
			next(file_in)
		for line in file_in:
			chain_ID.append(line[0])
			bturn_sequence.append(line[6]+line[13]+line[20]+line[27]) 
			#pos_seq+= line[2:5]+line[9:12]+line[16:19]+line[23:26]
			bturn_type.append(line[29:33])
			for i in(line[1:5],line[8:12],line[15:19],line[22:26]):
				aa_pos.append(i)
				

	#writting sequences in format file
	bturn_parsed_file = bturn_file[0:4]+"Output_bturn.fasta"
	with open(bturn_parsed_file,"w+") as file_out:
		file_out.write("> {} PARSED\n\n".format(bturn_file))
		file_out.write("            Postion                  Chain     Sequence\n")
		lign_number = len(bturn_sequence)//4
		j = 0
		for i in range(lign_number):
			file_out.write(" {}         {}          {}\n".format(aa_pos[j:j+4],chain_ID[i],bturn_sequence[j:j+4]))
			j += 4

	return bturn_sequence, aa_pos , bturn_type,chain_ID


"""def main():
		bturn_file = "1oip.bturns"
		btun_sequence, aa_pos,bturn_type = parse_promotif_bturn(bturn_file)
		print(btun_sequence)
		print(aa_pos)
		print(bturn_type)
	
main()
exit(0)"""











