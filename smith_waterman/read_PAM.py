import re

def read_matrix(filename):

	scoring_dict = {}
	flag = False
	with open(filename) as fh:
		count = 0
		for line in fh:
			line = line.strip()
			if line.startswith("#"): continue
			elif line.startswith("A"):
				amino_acids = line.split("  ")
				flag = True
				continue

			if flag:
				values = []
				m = re.split(' ',line) #could not figure out regex for 1 and 2 spaces
				for x in m:
					try: values.append(int(x))
					except: pass

				for index,value in enumerate(values):
					scoring_dict[(amino_acids[count],amino_acids[index])] = int(value)
				count += 1

	return scoring_dict

def read_optimized_matrix(filename):

	scoring_dict = {}
	flag = False
	AA = "ARNDCQEGHILKMFPSTWYV"
	with open(filename) as fh:
		for line in fh:
			line = line.strip()

			line = line.split(" ")
			aa = line[0]
			values = [float(x) for x in line[1::]]

			for index,value in enumerate(values):
				#print(index,value)
				scoring_dict[(aa,AA[index])] = float(value)
				scoring_dict[(AA[index],aa)] = float(value)
				scoring_dict[("X",aa)] = -1
				scoring_dict[(aa,"X")] = -1
				scoring_dict[("X","X")] = -1

				scoring_dict[("B",aa)] = -1
				scoring_dict[(aa,"B")] = -1
				scoring_dict[("B","B")] = -1

				scoring_dict[("Z",aa)] = -1
				scoring_dict[(aa,"Z")] = -1
				scoring_dict[("Z","Z")] = -1

	return scoring_dict

if __name__ == "__main__":
	print(read_matrix("../BLOSUM62"))




