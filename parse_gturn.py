def parse_promotif_gturn(gturn_file):
	gturn_sequence = []
	aa_pos = []
	chain_ID = []
	gturn_type = []
	#pos_seq = ""
	with open(gturn_file, 'r') as file_in:
		for i in range(2):
			next(file_in)

		for line in file_in:
			chain_ID .append(line[0])
			gturn_sequence.append(line[6]+line[13]+line[20]) 
			gturn_type.append(line[21:29])
			#pos_seq+= line[2:5]+line[9:12]+line[16:19]+line[23:26]
			for i in(line[2:5],line[9:12],line[16:19]):
				aa_pos.append(i)


	#writting sequences in format file
	gturn_parsed_file = gturn_file[0:4]+"Output_gturn.txt"
	with open(gturn_parsed_file,"w+") as file_out:
		file_out.write("> {} PARSED\n\n".format(gturn_file))
		file_out.write("            Position          Chain     Sequence      Type\n")
		lign_number = len(gturn_sequence)//3
		j = 0
		for i in range(lign_number):
			file_out.write("   {}        {}         {}      {}\n".format(aa_pos[j:j+3],chain_ID[i],gturn_sequence[j:j+3],gturn_type[i]))
			j += 3
	return gturn_sequence, aa_pos, gturn_type,chain_ID


"""def main():
		gturn_file = "1oip.gturns"
		gturn_sequence, aa_pos ,gturn_type,chain_ID= parse_promotif_gturn(gturn_file)
		print(gturn_sequence)
		print(aa_pos)
	
	
main()
exit(0)"""











