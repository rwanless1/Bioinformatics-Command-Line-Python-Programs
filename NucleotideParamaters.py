import sys





class FastAreader:
    '''
    Define objects to read FastA files.

    instantiation:
    thisReader = FastAreader ('testTiny.fa')
    usage:
    for head, seq in thisReader.readFasta():
        print (head,seq)
    '''

    def __init__(self, fname=None):
        '''contructor: saves attribute fname '''
        self.fname = fname

    def doOpen(self):
        ''' Handle file opens, allowing STDIN.'''
        if self.fname is None:
            return sys.stdin
        else:
            return open(self.fname)

    def readFasta(self):
        ''' Read an entire FastA record and return the sequence header/sequence'''
        header = ''
        sequence = ''

        with self.doOpen() as fileH:

            header = ''
            sequence = ''

            # skip to first fasta header
            line = fileH.readline()
            while not line.startswith('>'):
                line = fileH.readline()
            header = line[1:].rstrip()

            for line in fileH:
                if line.startswith('>'):
                    yield header, sequence
                    header = line[1:].rstrip()
                    sequence = ''
                else:
                    sequence += ''.join(line.rstrip().split()).upper()

        yield header, sequence


class NucParams:
    rnaCodonTable = {
        # RNA codon table
        # U
        'UUU': 'F', 'UCU': 'S', 'UAU': 'Y', 'UGU': 'C',  # UxU
        'UUC': 'F', 'UCC': 'S', 'UAC': 'Y', 'UGC': 'C',  # UxC
        'UUA': 'L', 'UCA': 'S', 'UAA': '-', 'UGA': '-',  # UxA
        'UUG': 'L', 'UCG': 'S', 'UAG': '-', 'UGG': 'W',  # UxG
        # C
        'CUU': 'L', 'CCU': 'P', 'CAU': 'H', 'CGU': 'R',  # CxU
        'CUC': 'L', 'CCC': 'P', 'CAC': 'H', 'CGC': 'R',  # CxC
        'CUA': 'L', 'CCA': 'P', 'CAA': 'Q', 'CGA': 'R',  # CxA
        'CUG': 'L', 'CCG': 'P', 'CAG': 'Q', 'CGG': 'R',  # CxG
        # A
        'AUU': 'I', 'ACU': 'T', 'AAU': 'N', 'AGU': 'S',  # AxU
        'AUC': 'I', 'ACC': 'T', 'AAC': 'N', 'AGC': 'S',  # AxC
        'AUA': 'I', 'ACA': 'T', 'AAA': 'K', 'AGA': 'R',  # AxA
        'AUG': 'M', 'ACG': 'T', 'AAG': 'K', 'AGG': 'R',  # AxG
        # G
        'GUU': 'V', 'GCU': 'A', 'GAU': 'D', 'GGU': 'G',  # GxU
        'GUC': 'V', 'GCC': 'A', 'GAC': 'D', 'GGC': 'G',  # GxC
        'GUA': 'V', 'GCA': 'A', 'GAA': 'E', 'GGA': 'G',  # GxA
        'GUG': 'V', 'GCG': 'A', 'GAG': 'E', 'GGG': 'G'  # GxG
    }
    aaComp = {
        'A': 0, 'G': 0, 'M': 0, 'S': 0, 'C': 0,
        'H': 0, 'N': 0, 'T': 0, 'D': 0, 'I': 0,
        'P': 0, 'V': 0, 'E': 0, 'K': 0, 'Q': 0,
        'W': 0, 'F': 0, 'L': 0, 'R': 0, 'Y': 0
    }
    nucComp = {}
    dnaCodonTable = {key.replace('U', 'T'): value for key, value in rnaCodonTable.items()}


#defines all the dictionaries which will be used by the methods in the program and the variables which count the values in the dictionary
    def __init__(self, inString=''):
        self.codonComp = {codon:0 for codon in NucParams.rnaCodonTable}
        self.aaComp = {aa:0 for aa in NucParams.rnaCodonTable.values()}
        self.nucComp = {nuc:0 for nuc in 'AGCTNU'}
        self.DNAcodonComp = {dna:0 for dna in NucParams.dnaCodonTable}
        self.nucCount = 0
        self.GCpercent = 0
        self.MbLength = 0


# gets the sequnce input and changes it from lowercale to uppercase and converts 'T' to 'U'
    def addSequence(self, inSeq):
        index = 0
        upper = inSeq.upper()
        rnaseq = upper.replace('T', 'U')
        codon = rnaseq[index: index + 3]
#counts the number of nuclotided and adds it to nuComp dictionary
        for codon in upper:
            self.nucComp[codon] += 1
# sets the codon length at three and counts every three base pairs
        seqlen = len(upper)
        for index in range(0, seqlen, 3):
            if 'T'in upper:
                rnaseq = upper.replace('T', 'U')
            #sets codon length at 3 and then adds a new codon for every 3rd string input
            codon = rnaseq[index: index + 3]
            #adds nucleotide and codon count with input string
            # matches nucComp values to corresponding amino acid values
            # and counts for amino acids
            if codon in self.rnaCodonTable:
                self.codonComp[codon] += 1
                aa = self.rnaCodonTable.get(codon)
                self.aaComp[aa] += 1
            else:
                pass
        self.nucCounter = sum(self.nucComp.values())



        #calculates gc number by sum of G and C / (total length of sequence)
        gValue = int(self.nucComp.get('G'))
        cValue = int(self.nucComp.get('C'))

        gcSum = gValue + cValue
        self.GCpercent = ((gcSum)/(self.nucCounter) * 100)


        #calculates how many megabases in the sequence by calculated total length of sequnce by nuccounter by 1million

        self.MbLength = (self.nucCounter / 1000000)

        return

    def aaComposition(self, inString):
        # returns the dictionary with aaComp values inputed
        return self.aaComp

    def nucComposition(self):
        #returns the nunComp dictionary with the calculated values

        return self.nucComp

    def codonComposition(self):
        return self.codonComp

    def nucCount(self):
        self.nucCounter = 0
        # calculates the total of all the nucleotides in nucComp Dictionary
        self.nucCounter = sum(self.nucComp.values())

        return self.nucCounter





def main(fileName=None):
    myReader = FastAreader(fileName)
    myNuc = NucParams()


    for head, seq in myReader.readFasta():
        myNuc.addSequence(seq)

        # sort codons in alpha order, by Amino Acid
    print ("Sequence length =", f'{myNuc.MbLength: .2f}', "Mb" )

    print("GC content: ", f'{myNuc.GCpercent: .1f}', "%")

        # calculate relative codon usage for each codon and print
    for codon in myNuc.codonComp.items():
        print('{:s} : {:s} {:5.1f} ({:6d})'.format(codon[0], myNuc.rnaCodonTable[codon[0]], (codon[1] / sum (myNuc.codonComp.values())) * 100, codon[1]))

if __name__ == "__main__":
     main()  # make sure to change this in order to use stdin
