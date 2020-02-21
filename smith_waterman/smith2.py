from PAM import PAM
from read_PAM import read_matrix


class local_alignment (object):
    """  
        local_alignment obviously my doods
    """

    def __init__(self, query, target, gap_penalty=-5, gap_extension=-2, score_dict=None):


        self.gap_penalty = gap_penalty
        self.gap_extension = gap_extension

        
        self.query = query    
        self.target = target    

        self.querylen = len(query)
        self.targetlen = len(target)
        self.score_dict = score_dict

        self.table = [[0 for i in range(self.targetlen)] for j in range(self.querylen)]

        self.traceback = [[(0, 0) for i in range(self.targetlen)] for j in range(self.querylen)]


        return

    def score(self):    

        highscore = high_i = high_j = 0       # highest scores encountered so far in the matrix

        best_q_alignment = []  
        best_t_alignment = []  

        for i in range(self.querylen):  # score the left edge


            self.table[i][0] = self.score_dict[(self.query[i],self.target[0])]

        for j in range(self.targetlen):  # score the top edge

            self.table[0][j] = self.score_dict[(self.target[j],self.query[0])]


        # Nested loop through the remainder of the matrix, scoring as we go
        prev_t_gap = False
        prev_q_gap = False

        for i in range(1, self.querylen):      # start these iterations from 1, not 0, as we have already done edges

            for j in range(1, self.targetlen):
                #i = 2 and j = 2
                queryword = self.query[i: i + 1]  # An array slice is perhaps more natural in python than a substring
                targetword = self.target[j: j + 1]

                increment = self.score_dict[(queryword,targetword)]

                matchscore = self.table[i - 1][j - 1] + increment   # increment will contain either a positive reward
                                                                # or a negative penalty depending on whether we matched
                if prev_t_gap:
                    target_gap_score = self.table[i][j - 1] + self.gap_extension    # scores associated with gapping
                elif prev_q_gap:
                    query_gap_score = self.table[i - 1][j] + self.gap_extension
                else:
                    target_gap_score = self.table[i][j - 1] + self.gap_penalty    # scores associated with gapping
                    query_gap_score = self.table[i - 1][j] + self.gap_penalty     # in either the target or query

                best_score = max(
                    (0, (0, 0)),                        # 0 score will never have a traceback
                    (matchscore, (1, 1)),               # A match corresponds to a -1,-1 traceback
                    (target_gap_score, (0, 1)),         # A target gap corresponds to a 0, -1 traceback
                    (query_gap_score, (1, 0))           # A query gap corresponds to a -1, 0 traceback

                )

                if target_gap_score == best_score[0]:
                    prev_t_gap = True
                elif query_gap_score == best_score[0]:
                    prev_q_gap = True
                else:
                    prev_t_gap = False
                    prev_q_gap = False
                    

                self.table[i][j] = best_score[0]    # The first element in the tuple is the actual score to be recorded
                self.traceback[i][j] = best_score[1]    # The traceback offsets associated with the score are in a tuple


                if self.table[i][j] > highscore:    # This represents the "high road" approach.
                                                    # "low road" would be >=

                    highscore = self.table[i][j]    # record the new high score
                    high_i = i                      
                    high_j = j


        i = high_i          
        j = high_j

        while self.table[i][j] and i > -1 and j > -1:

            i_offset, j_offset = self.traceback[i][j]       # unpack the offset tuples stored in the traceback table

            if i_offset:
                best_q_alignment.append(self.query[i])
            else:
                best_q_alignment.append('-')                # if the value is a zero, we are gapping!

            if j_offset:

                best_t_alignment.append(self.target[j])

            else:
                best_t_alignment.append('-')                # if the value is a zero, we are gapping, now the other way

            i -= i_offset
            j -= j_offset
            # print(i,j)

            if i == 0 or j == 0:
                break

        best_q_alignment.reverse()  # flip 'em both once we are done, since we built them "end-to-beginning"
        best_t_alignment.reverse()

        out_string = '\nAlignment score ' + str(highscore) + ' and is:\n\nTarget:\t' + \
            str(j + 2) + '\t' + ''.join(best_t_alignment) + '\n\t\t\t'

        for k in range(len(best_t_alignment)):     # t and q alignments should be the same length!

            if best_t_alignment[k] == best_q_alignment[k]:

                out_string += ''    # Only put a bar if the two characters are identical at this position

            else:

                out_string += ' '    # otherwise just insert a space

        out_string += '\nQuery:\t' + str(i + 2) + '\t' + ''.join(best_q_alignment) + '\n'


        #print(highscore)
        #return out_string
        return highscore

    def __str__(self):                         
        lineout = 'Scoring table:\n\t' + '\t'.join(self.target) + '\n'
        # # The above is just a fancy looking way to break the target string into tab-delimited individual characters

        for i in range(self.querylen):
            lineout += self.query[i] + "\t"
            for j in range(self.targetlen):

                lineout += str(self.table[i][j]) + "\t"

            lineout += '\n'
        lineout += '\n\nTraceback table:\n\t' + '\t'.join(self.target) + '\n'

        for i in range(self.querylen):

            lineout += self.query[i] + "\t"

            for j in range(self.targetlen):

                lineout += ''.join([str(k) for k in self.traceback[i][j]]) + "\t"  # prettying up the traceback tuples

            lineout += '\n'

        return lineout


if __name__ == "__main__":
    # B = PAM(N=120)
    # PAM1 = B.Build_PAMN()

    BLOSUM62 = read_matrix("../BLOSUM62")

    seq1 = "MDSVCPQGKYIHPQNNSICCTKCHKGTYLYNDCPGPGQDTDCRECESGSFTASENHLRHC"

    seq2 = "LSCSKCRKEMGQVEISSCTVDRDTVCGCRKNQYRHYWSENLFQC"
    
    A = local_alignment(seq1, seq2,-8,-3,BLOSUM62) #pass dictionary

    print(A.score())

