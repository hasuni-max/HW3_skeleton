import numpy as np
import tqdm as tq
import os
import matplotlib.pyplot as plt
from smith_waterman import training_files, smith3, read_PAM, ROC, PAM, optimization

def test_roc():

	pos_matches = []
	all_files = []
	for file in training_files.grab_pairs("Pospairs.txt"):
		pos_matches.append(file)
		all_files.append(file[0])
		all_files.append(file[1])

	neg_matches = []
	for file in training_files.grab_pairs("Negpairs.txt"):
		neg_matches.append(file)
		all_files.append(file[0])
		all_files.append(file[1])

	all_files = set(all_files)

	sequences = {file:training_files.parse_fasta(file) for file in all_files}

	filepath = os.path.join("LAK_optimized")

	ah = read_PAM.read_optimized_matrix(filepath)
	matrices = [ah]

	thresholds = [x for x in np.arange(20, 21, 0.2)]
	gap_p = -11.6
	gap_e = -5.7
	r = ROC.roc()

	tps = []
	fps = []

	for index,matrix in enumerate(matrices):
		print("Testing LAK")

		for threshold in tq.tqdm(thresholds):##loop over 8 different thresholds
			## print(threshold)
			tp,fp = optimization.calculate_tp_fp(pos_matches,neg_matches,sequences,
				threshold,gap_p,gap_e,matrix, normalize = False)
			tps.append(tp)
			fps.append(fp)

	assert len(tps) == len(fps)


def test_smithwaterman():

	BLOSUM50 = read_PAM.read_matrix("BLOSUM50")
    
	seq1 = "MDSVCPQGKYIHPQNNSICCTKCHKGTYLYNDCPGPGQDTDCRECESGSFTASENHLRHC"

	seq2 = "LSCSKCRKEMGQVEISSCTVDRDTVCGCRKNQYRHYWSENLFQC"
	
	A = smith3.local_alignment(seq1, seq2,-11,-1,BLOSUM50)
	assert A.score() > 0

def test_perfect_alignment():



	BLOSUM50 = read_PAM.read_matrix("BLOSUM50")
    
	seq1 = "LSCSKCRKEMGQVEISS"

	seq2 = "LSCSKCRKEMGQVEISS"
	
	A = smith3.local_alignment(seq1, seq2,-11,-1,BLOSUM50)
	assert A.score() > 0

def test_scoring():
    return None
