import tqdm as tq
import numpy as np
from .training_files import grab_pairs, parse_fasta
from .smith2 import local_alignment 
from .read_PAM import read_optimized_matrix
from .PAM import PAM
from .ROC import roc


def run_local_alignment(seq1, seq2, gap_penalty, gap_extension,matrix):
	A = local_alignment(seq1.upper(), seq2.upper(),gap_penalty,gap_extension,matrix) #pass dictionary
	return A.score()

def calculate_tp_fp(pos_matches,neg_matches,sequences,threshold,gap_p,gap_e,matrix,normalize=False):
	"""
		Calculate true positives and false positives given positive and negative
		matches, their corresponding sequences (sequences), a gap opening and a 
		gap extension penalty and finally a scoring matrix
	"""
	true_pos = 0
	false_neg = 0 
	true_neg = 0
	false_pos = 0 


	for pos in pos_matches:
		score = run_local_alignment(sequences[pos[1]], sequences[pos[0]],gap_p,gap_e,matrix)

		if normalize:
			score = score/min(len(sequences[pos[1]]),len(sequences[pos[0]]))

		if score >= threshold:
			true_pos += 1
		else:
			false_neg += 1	
			
	for neg in neg_matches:
		score = run_local_alignment(sequences[neg[1]], sequences[neg[0]],gap_p,gap_e,matrix)

		if normalize:
			score = score/min(len(sequences[neg[1]]),len(sequences[neg[0]]))

		if score < threshold:
			true_neg += 1
		else:
			false_pos += 1	

	tp = true_pos/(true_pos+false_neg)
	fp = false_pos/(true_neg+false_pos)	


	return tp,fp



def adjust_PAM(PAM):

	AA = "ARNDCQEGHILKMFPSTWYV"
	for aa in AA:
		PAM[("X",aa)] = -1
		PAM[(aa,"X")] = -1
		PAM[("X","X")] = -1
		PAM[("B",aa)] = -1
		PAM[(aa,"B")] = -1
		PAM[("B","B")] = -1
		PAM[("Z",aa)] = -1
		PAM[(aa,"Z")] = -1
		PAM[("Z","Z")] = -1

	return PAM

	

if __name__ == "__main__":


	pos_matches = []
	all_files = []
	for file in grab_pairs("../Pospairs.txt"):
		pos_matches.append(file)
		all_files.append(file[0])
		all_files.append(file[1])

	neg_matches = []
	for file in grab_pairs("../Negpairs.txt"):
		neg_matches.append(file)
		all_files.append(file[0])
		all_files.append(file[1])

	all_files = set(all_files)

	sequences = {file:parse_fasta("../"+file) for file in all_files}



	############################ nPAM matrices ###############################
	matrices = []
	labels = []
	for n in range(20,320,60):
		labels.append(str(n))
		A = PAM(N=n)
		new_PAM = A.Build_PAMN()
		matrices.append(adjust_PAM(new_PAM))
		##break
	


	## Make ROC curves. Run once with normalize False and True
	thresholds = [x for x in range(0,300,10)]
	gap_p = -8
	gap_e = -3
	r = roc()

	for index,matrix in enumerate(matrices):
		print("Testing",labels[index])

		for threshold in tq.tqdm(thresholds):#loop over 8 different thresholds
			## print(threshold)
			tp,fp = calculate_tp_fp(pos_matches,neg_matches,sequences,
				threshold,gap_p,gap_e,matrix, normalize = False)
			## print(tp,fp)

			r.add_rates(tp,fp)
		r.plot_ROC(lab=labels[index])
		r.new_curve()

	r.save_plot("nPAM_ROC")




	########################################## LAK #########################################

	ah = read_optimized_matrix("../LAK_optimized")
	matrices = [ah]

	thresholds = [x for x in np.arange(20, 60, 0.2)]
	gap_p = -11.6
	gap_e = -5.7
	r = roc()


	for index,matrix in enumerate(matrices):
		print("Testing LAK")

		for threshold in tq.tqdm(thresholds):##loop over 8 different thresholds
			## print(threshold)
			tp,fp = calculate_tp_fp(pos_matches,neg_matches,sequences,
				threshold,gap_p,gap_e,matrix, normalize = False)
			## print(tp,fp)

			r.add_rates(tp,fp)
		r.plot_ROC(lab="LAK_optimized")
		r.new_curve()

	r.save_plot("LAK_optimized")




	






