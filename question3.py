'''
Counts nucleotides in DNA sequence
'''
class dnaString (str):
    def length (self):
        return (len(self))

    def countNucleotideA (self):
        num_A = self.count('A')
        return num_A
    def countNucleotideT (self):
        num_T = self.count('T')
        return num_T
    def countNucleotideG (self):
        num_G = self.count('G')
        return num_G
    def countNucleotideC (self):
        num_C = self.count('C')
        return num_C

dna = input("Enter a dna sequence: ")
upperDNA = dna.upper()
coolString = dnaString(upperDNA)

print ("number As = {} number Cs = {} number Gs = {} number Ts = {}".format(
    coolString.countNucleotideA(),
    coolString.countNucleotideC(),coolString.countNucleotideG(), coolString.countNucleotideT() ) )