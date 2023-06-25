'''
Counts A and T in DNA sequence
'''
class dnaString (str):
    def length (self):
        return (len(self))

    def getAT (self):
        num_A = self.count('A')
        num_T = self.count('T')
        return ((num_A + num_T)/ self.length() ) # self.length was misspelled

dna = input("Enter a dna sequence: ")
upperDNA = dna.upper()
coolString = dnaString(upperDNA)

print ('AT content = {0:0.3f}''.format(coolString.getAT()) ) #the {0:0.3f} was used to have 3 decimal places
