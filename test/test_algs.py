import numpy as np
import tqdm as tq
import os
import matplotlib.pyplot as plt
from smith_waterman import training_files, smith2, read_PAM, ROC, PAM,optimization

def test_roc():

	filepath = os.path.join("HW3_skeleton", "LAK_optimized")

	ah = read_PAM.read_optimized_matrix(filepath)
	matrices = [ah]

	thresholds = [x for x in np.arange(20, 60, 0.2)]
	gap_p = -11.6
	gap_e = -5.7
	r = roc()

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
			## print(tp,fp)

			r.add_rates(tp,fp)
		r.plot_ROC(lab="Test")
		r.new_curve()


	assert len(tp) == len(fp)


def test_smithwaterman():

	BLOSUM50 = read_PAM.read_matrix("HW3_skeleton/BLOSUM50")
    
	seq1 = "MDSVCPQGKYIHPQNNSICCTKCHKGTYLYNDCPGPGQDTDCRECESGSFTASENHLRHC"

	seq2 = "LSCSKCRKEMGQVEISSCTVDRDTVCGCRKNQYRHYWSENLFQC"
	
	A = smith_waterman.local_alignment(seq1, seq2,-11,-1,BLOSUM50)
	assert A.score() > 0

def test_scoring():
    return None
