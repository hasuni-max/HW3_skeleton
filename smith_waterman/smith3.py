from .PAM import PAM
from .read_PAM import read_matrix


class local_alignment (object):
    """  
        An object that performs the standard smithwaterman alignment
    """

    def __init__(self, query, target, gap_p=-5, gap_e=-2, score_dict=None):


        self.gap_p = gap_p 
        self.gap_e = gap_e

        
        self.query = query    
        self.target = target    

        self.querylen = len(query)
        self.targetlen = len(target)
        self.score_dict = score_dict #This contains the substitution matrix that we pass when we call the class


        #initialize an empty table with zeros. The table has the dimension querylen X targetlen
        #This will contain all the scores of the matrix
        #Will be filled when object method score() is called
        self.table = [[0 for i in range(self.targetlen)] for j in range(self.querylen)]

        #initialize an empty table with zeros. The table has the dimension querylen X targetlen
        #This will contain all the traceback informaton
        #Will be filled when object method score() is called
        self.traceback = [[(0, 0) for i in range(self.targetlen)] for j in range(self.querylen)]


        return

    def score(self):    

        #These variables keep track of the highest alignment score and its position in the
        #scoring table and the tracback table. These scors will be updated as the we traverse both sequences 
        highscore = high_i = high_j = 0       


        #Once both sequenes have been traversed and all scoeres have been calculated. 
        #we loop back through and add amino acids back
        best_q_alignment = []  
        best_t_alignment = []  


         # score the left edge. It is easier to do this now. We do the same thing with the right edge
        for i in range(self.querylen): 
            self.table[i][0] = self.score_dict[(self.query[i],self.target[0])]

        #score the right edge
        for j in range(self.targetlen):  # score the top edge
            self.table[0][j] = self.score_dict[(self.target[j],self.query[0])]


        #This dictionary keeps track of gaps.
        gapsies = {}

        #This is the main body of the code. We simply loop through the strings and score as we go, 
        #keeping track of the highest score as we go along.
        for i in range(1, self.querylen):      # start these iterations from 1, not 0, as we have already done edges

            for j in range(1, self.targetlen):

                #Current amino acid from query and target
                queryword = self.query[i: i + 1]  #string slice
                targetword = self.target[j: j + 1]

                increment = self.score_dict[(queryword,targetword)] #grab score from substitution matrix

                #add the increment to the value from the diagonal 
                matchscore = self.table[i - 1][j - 1] + increment   # The increment from above either has a 
                #positive or negative reward
                         

                #REMOVE FLAG and added this dictionary that keeps track gaps
                if gapsies.get((i,j-1),0) == 1:
                    target_gap_score = self.table[i][j - 1] + self.gap_e 
                else:
                    target_gap_score = self.table[i][j - 1] + self.gap_p

                if gapsies.get((i-1,j),0) == -1:
                    query_gap_score = self.table[i - 1][j] + self.gap_e
                else:
                    query_gap_score = self.table[i - 1][j] + self.gap_p
                

                #So this is a nifty Python trick. Basically checks the first index in each tuple and takes the max
                #The value returned is the tuple associated with the max first index
                #The nested loops contain values for the traceback table
                best_score = max(
                    (0, (0, 0)),                        # 0 score will never have a traceback
                    (matchscore, (1, 1)),               # A match corresponds to a -1,-1 traceback
                    (target_gap_score, (0, 1)),         # A target gap corresponds to a 0, -1 traceback
                    (query_gap_score, (1, 0))           # A query gap corresponds to a -1, 0 traceback

                )

                #Update the gap flag. If the best score from previous step was a gap score then set flag to true
                #else we set flag to false
                if target_gap_score == best_score[0]:
                    gapsies[(i,j)] = 1
                elif query_gap_score == best_score[0]:
                    gapsies[(i,j)] = -1
                
                #Update the scoring table and the traceback table 
                self.table[i][j] = best_score[0]    # The first element in the tuple is the best score and added to the table
                self.traceback[i][j] = best_score[1]   #Second element is the tuple with the traceback back values 

                #If the current position in the table is greater than the highscore
                #If so update highscore and highscore positions
                if self.table[i][j] > highscore:   
                                                    # "low road" would be >=

                    highscore = self.table[i][j]    # record the new high score
                    high_i = i                      
                    high_j = j


        #Loop has finished. Set i and j to the indices associated with the highest scores
        i = high_i          
        j = high_j

        #This loop grabs the sequence associated with the best alignment
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

        #Reverse the sequences 
        best_q_alignment.reverse()  
        best_t_alignment.reverse()


        #The rest of the method creates a nicely formatteed alignment output string. This can be printed
        #for an obect n by typing n.out_string
        self.out_string = '\nAlignment score ' + str(highscore) + ' and is:\n\nTarget:\t' + \
            str(j + 2) + '\t' + ''.join(best_t_alignment) + '\n\t\t\t'

        for k in range(len(best_t_alignment)):     # t and q alignments should be the same length!

            if best_t_alignment[k] == best_q_alignment[k]:

               self.out_string += ''    # Only put a bar if the two characters are identical at this position

            else:

                self.out_string += ' '    # otherwise just insert a space

        self.out_string += '\nQuery:\t' + str(i + 2) + '\t' + ''.join(best_q_alignment) + '\n'


        #print(highscore)
        #return out_string
        return highscore

    def __str__(self):

        """
            This class has a string attribute that will print the scoring table and the traceback table
        """                     
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

    seq1 = "LSCSKCRKEMGQVEISSCTVDRDTVCGCRKNQYRHYWSENLFQC"

    seq2 = "LSCSKCRKEMGQVEISSCTVDRDTVCGCRKNQYRHYWSENLFQC"
    
    A = local_alignment(seq1, seq2,-11,-1,BLOSUM62) #pass dictionary

    print("Alignment score: ",A.score())



