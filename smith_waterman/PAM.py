from math import log10


class PAM(object):

    """
        This class calculates multiples of a PAM matrix starting with a PAM1. 
    """

    def __init__(self, N=250, filename='PAM1.txt', aa_frequencies=None):

        self.PAM = {}
       
        self.N = N                    # The evolutionary distance of the PAM matrix we desire
        self.filename = filename      

        #These are background AA frequencies found on the internet. I am not sure if these values are accurate 
        if aa_frequencies is None:

            self.aa_frequencies = {"G": 0.089, "R": 0.041, "A": 0.087, "N": 0.040, "L": 0.085, "F": 0.040, "K": 0.081,
                                   "Q": 0.038, "S": 0.070, "I": 0.037, "V": 0.065, "H": 0.034, "T": 0.058, "C": 0.033,
                                   "P": 0.051, "Y": 0.030, "E": 0.050, "M": 0.015, "D": 0.047, "W": 0.010}
        else:

            self.aa_frequencies = aa_frequencies  # OK, so the normalized amino acid frequencies are a gimme too!

    def Build_PAMN(self):

        PAM1, alphabet = self.read_PAM1()       # Read in the PAM1 file

        size = len(alphabet)  # Record how many aa were in the input file
        # Some PAM table include ambiguity characters, so we cannot assume 20

        PAM1 = [[element / 10000 for element in row] for row in PAM1]


        PAMN = list(PAM1)    


        #Keep multiplying matrix until N - 1.... we subtract one because we start with PAM1
        for i in range(self.N - 1):  

            PAMN = self.__MatrixMultiply(PAMN, PAM1)

        for i in range(size):  # Convert to a log odds formulation.
            for j in range(size):
                PAMN[i][j] = log10(PAMN[i][j] / self.aa_frequencies[alphabet[i]])

        for i in range(size):  
            for j in range(size):
                self.PAM[alphabet[i], alphabet[j]] = int(round(10 * ((PAMN[i][j] + PAMN[j][i]) / 2)))

        return self.PAM

    def read_PAM1(self):

        """Read in the raw PAM1 file into a two-dimensional list"""

        PAM1 = []       

        with open(self.filename, 'r') as pam1_file:

            line = pam1_file.readline()  
            line = line.strip()          

            alphabet = line.split(' ')

            del alphabet[0]   
            line = pam1_file.readline()    

            while line:

                line = line.strip()                             
                line = line.split(' ')                         
                del line[0]                                     
                line = [int(element) for element in line]       
                PAM1.append(line)
                line = pam1_file.readline()

        return PAM1, alphabet

    def PAM_dump(self):

        """Print out the final scoring dict we have created. """

        ordered_symbols = self.aa_frequencies.keys()       
        ordered_symbols.sort()

        for aa in ordered_symbols:
            print('\t', aa,)
        print('\n')
        for aa1 in ordered_symbols:
            print(aa1, '\t',)
            for aa2 in ordered_symbols:
                print(self.PAM[aa1, aa2], '\t',)
            print('\n')

    def __MatrixMultiply(self, arrayA, arrayB):

        """Multiply two matrices, A and B, and return matrix C. May dot products haunt your dreams!"""

        rowsA = len(arrayA)         # Rows of A
        colsB = len(arrayB[0])      # Columns of B

        if rowsA != colsB:

            raise ValueError("math domain error")

        new_array = [[0 for i in range(rowsA)] for j in range(colsB)]  # initialize an empty array of correct size

        rowcols = len(arrayA[0])                         # Number of columns in A and rows in B

        for i in range(rowsA):                          # over rows of A
                for j in range(colsB):                  # over columns of B.
                        for k in range(rowcols):        # over columns of A and rows of B

                                new_array[i][j] += arrayA[i][k] * arrayB[k][j]

        return new_array


def main():

    A = PAM(N=200)
    print(A.Build_PAMN())


    #A.PAM_dump()

if __name__ == '__main__':

    main()