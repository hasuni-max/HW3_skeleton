import sys
import matplotlib.pyplot as plt
import numpy as np
from .smith3 import local_alignment
from .read_PAM import read_matrix


BLOSUM62 = read_matrix("BLOSUM62")

seq1 = "MDSVCPQGKYIHPQNNSICCTKCHKGTYLYNDCPGPGQDTDCRECESGSFTASENHLRHC"

seq2 = "LSCSKCRKEMGQVEISSCTVDRDTVCGCRKNQYRHYWSENLFQC"

A = local_alignment(seq1, seq2,-8,-3,BLOSUM62) #pass dictionary

print("Alignment score: ",A.score())

print(A.out_string)