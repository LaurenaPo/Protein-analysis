def parse_promotif_bbulge(bbulge_file):
	x_seq = []
	first_seq = []
	second_seq = []
	x_pos = []	
	first_pos = []
	second_pos = []
	chain_ID = []
	bbluge_type = []
	with open(bbulge_file, 'r') as file_in:
		#line = file_in.readline().split()
		#line = file_in.readline().strip()	
		for line in file_in:
			x_seq.append(line[11])
			first_seq.append(line[35])
			second_seq.append(line[59])
			x_pos.append(line[5:9])
			first_pos.append(line[29:34])
			second_pos.append(line[53:57])
			bbluge_type.append(line[0:4])


	#writting sequences in format file
	bbulge_parsed_file = bbulge_file[0:4]+"Output_bbulge.fasta"
	with open(bbulge_parsed_file,"w+") as file_out:
		file_out.write("> {} PARSED\n\n".format(bbulge_file))
		#Column 1: position
		#column 2 : sequence
		# column 3 : bbugle type
		file_out.write("    Position              Sequence       type\n")
		#x : residu on the normal strend
		#1 and 2 : residus on the bulged strands
		file_out.write("  x    1      2        x    1    2       \n\n")
		for i in range(len(x_seq)):
			file_out.write("{}  {} {}       {}    {}    {}      {}\n".format(x_pos[i],first_pos[i],second_pos[i],x_seq[i],first_seq[i],second_seq[i],bbluge_type[i]))
		

	return x_seq, x_pos, first_seq, first_pos, second_seq, second_pos, bbluge_type, chain_ID

"""def main():
		bbulge_file = "1oip.blg"
		x_seq, x_pos, first_seq, first_pos, second_seq, second_pos, bbluge_type, chain_ID = parse_promotif_bbulge(bbulge_file)
		print(x_seq)
		
main()
exit(0)"""











