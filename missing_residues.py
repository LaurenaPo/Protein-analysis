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




def compare_seq( pos_base, sequence_base, pos_test,sequence_test):
	#liste des positions ou la sequence n'est pas identique
	pos_difference=[]

	for i in range(len(pos_base)):
		if pos_base[i] == pos_test[i]:
			if sequence_base[i] != sequence_test[i]:
				pos_difference.append(i)
	return pos_difference 

"""def main():
    pos_base=['1','2','3','4']
    sequence_base='BBCD'
    pos_test=['1','2','3','4']
    sequence_test='ABCD'
    print(compare_seq(pos_base, sequence_base, pos_test,sequence_test))


main()"""
