import tqdm as tq
from training_files import grab_pairs, parse_fasta
from smith2 import local_alignment 
from read_PAM import read_matrix
from ROC import roc

# from smith_protein import local_alignment 


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
			#print(score)
		#print(score,threshold)
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

#Question 1
def find_best_gaps(pos_matches,neg_matches,sequences):

	"""
		Find the best gap penalty and gap extension combination.

		Uses BLOSUM50

	"""
	BLOSUM50 = read_matrix("../BLOSUM50")
	threshold = 115
	best_fp = 10
	for gap_p in tq.tqdm(range(1,21,1)): #try each gap penality from 1 to 20
		for gap_e in range(1,6,1): #try each gap extension from 1 to 5
			tp,fp = calculate_tp_fp(pos_matches,neg_matches,sequences,threshold,-gap_p,-gap_e,BLOSUM50)

			if fp < best_fp:
				best_fp = fp
				best_gap_p = gap_p
				best_gap_e = gap_e

	#print(best_gap_p,best_gap_e)
	return best_gap_p,best_gap_e




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



	BLOSUM50 = read_matrix("../BLOSUM50")

	# find_best_gaps(pos_matches,neg_matches,sequences)


	#Run the line below to find best gaps and neg matches. BEWARE WILL TAKE LONG
	#find_best_gaps(pos_matches,neg_matches,sequences)



	# Make ROC curves. Run once with normalize False and True
	matrices = ["BLOSUM50","BLOSUM62","PAM100","PAM250"]
	thresholds = [0,10,15,30,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,210,230,250,500]
	#thresholds = [.1,.2,.4,.6,.8,1,1.2,1.4,1.6,1.8,2,2.2,2.4]
	gap_p = -8
	gap_e = -3
	r = roc()
	for matrix in matrices:
		matrix_dict = read_matrix("../" + matrix)
		print("Testing ",matrix)
		for threshold in tq.tqdm(thresholds):#loop over 8 different thresholds
			
			tp,fp = calculate_tp_fp(pos_matches,neg_matches,sequences,
				threshold,gap_p,gap_e,matrix_dict, normalize = False)

			r.add_rates(tp,fp)
		r.plot_ROC(lab=matrix)
		r.new_curve()

	r.save_plot("final_ROC")






