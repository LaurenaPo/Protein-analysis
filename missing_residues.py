#recovery of the missing residues position
def missing_residues(pdb_file):
	missing_seq=[]
	with open(pdb_file, 'r') as file_in:
		for line in file_in:
			if line[0:27]=="REMARK 465   M RES C SSSEQI":
				test=line[0:10]
				line_i=file_in.readline()
				while(test == "REMARK 465"):
					missing_seq.append(line_i[22:26])
					line_i=file_in.readline()
					test=line_i[0:10]
	return missing_seq



