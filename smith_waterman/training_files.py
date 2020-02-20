import os

def grab_pairs(filename):

	with open(filename) as fh:
		for line in fh:
			line = line.strip().split()
			yield line[0], line[1]

def parse_fasta(filename):

	seq = ""
	with open(filename) as fh:
		for line in fh:
			if line.startswith(">"):
				continue
			seq += line.strip()

	return seq