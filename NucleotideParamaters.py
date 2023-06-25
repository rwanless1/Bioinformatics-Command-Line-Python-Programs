class ProteinParam:
    # These tables are for calculating:
    #     molecular weight (aa2mw), along with the mol. weight of H2O (mwH2O)
    #     absorbance at 280 nm (aa2abs280)
    #     pKa of positively charged Amino Acids (aa2chargePos)
    #     pKa of negatively charged Amino acids (aa2chargeNeg)
    #     and the constants aaNterm and aaCterm for pKa of the respective termini
    #  Feel free to move these to appropriate methods as you like

    # As written, these are accessed as class attributes, for example:
    # ProteinParam.aa2mw['A'] or ProteinParam.mwH2O

    aa2mw = {
        'A': 89.093, 'G': 75.067, 'M': 149.211, 'S': 105.093, 'C': 121.158,
        'H': 155.155, 'N': 132.118, 'T': 119.119, 'D': 133.103, 'I': 131.173,
        'P': 115.131, 'V': 117.146, 'E': 147.129, 'K': 146.188, 'Q': 146.145,
        'W': 204.225, 'F': 165.189, 'L': 131.173, 'R': 174.201, 'Y': 181.189
    }

    mwH2O = 18.015
    aa2abs280 = {'Y': 1490, 'W': 5500, 'C': 125}

    aa2chargePos = {'K': 10.5, 'R': 12.4, 'H': 6}
    aa2chargeNeg = {'D': 3.86, 'E': 4.25, 'C': 8.33, 'Y': 10}
    aaNterm = 9.69
    aaCterm = 2.34
    aaComp = {
            'A': 0, 'G': 0, 'M': 0, 'S': 0, 'C': 0,
            'H': 0, 'N': 0, 'T': 0, 'D': 0, 'I': 0,
            'P': 0, 'V': 0, 'E': 0, 'K': 0, 'Q': 0,
            'W': 0, 'F': 0, 'L': 0, 'R': 0, 'Y': 0
            }

    def __init__(self, protein):
        '''Build initial AA composition.'''

       
            # count symbols in protein
            # ignoring any bad characters
        for aa in protein.upper():
            if aa in self.aaComp:  # count valid AA
                self.aaComp[aa] += 1

    def aaCount(self):

        return sum (self.aaComp.values())
 # sum aa.

    # ceates dictionary named pI and goes in range of 1 to 1401, then divides by 100 and gets absolute value
    # then finds minimum of key value
    def pI(self):
        pI = {}
        for ph in range(0, 1401):
            ph = ph/100
            charge = self._charge_(ph)
            pI[str(ph)] = abs(charge)
        #print(min(pI, key=pI.get))
        return min(pI, key=pI.get)

    # takes the values of aaComp dictionary and divides them by total number of amino acids in table

    def aaComposition(self,):
        for key in self.aaComp.keys():
            percentage = (self.aaComp[key] / self.aaCount()*100) # then multiplies by 100 to get the percent and 0.2f for two decimal places
            print(key, '=','%.2f' % percentage,'%' )

    # this is an imput of the math equation for total charge of an amino acid
    # subtracts positive by negative to get toal charge
    def _charge_(self, pH):
        positiveCharge = 0
        negativeCharge = 0
        for aa in self.aa2chargePos: # this is positive charge
            aaCount = self.aaComp[aa]
            positiveCharge += aaCount * ((10**self.aa2chargePos[aa]) / (10**self.aa2chargePos[aa] + 10**pH))
        positiveCharge += (10**self.aaNterm) / (10**self.aaNterm + 10**pH)

        for aa in self.aa2chargeNeg: # this is negative charge
            aaCount = self.aaComp[aa]
            negativeCharge += aaCount * ((10**pH) / (10**self.aa2chargeNeg[aa] + 10**pH))
        negativeCharge += (10**pH) / (10**self.aaCterm + 10**pH)


        return positiveCharge - negativeCharge

    def molarExtinction(self):
        # multipies values of aa comp dictionary and aa2abs280 value with corresponsind amino acid
        aa2abs280 = {'Y': 1490, 'W': 5500, 'C': 125}
        yvalue = self.aaComp['Y'] * self.aa2abs280['Y']
        wValue = self.aaComp['W'] * self.aa2abs280['W']
        cValue = self.aaComp['C'] * self.aa2abs280['C']
        totalV = yvalue+wValue+cValue # total of w y and c value
        return totalV


    def massExtinction(self):
        myMW = self.molecularWeight()
        return self.molarExtinction() / myMW if myMW else 0.0

    def molecularWeight(self):
        # add up molecular weigth of all amino acids times the number of them

        mwH20 = 18.015
        molecularweight = mwH20
        for aa, mw in self.aa2mw.items():
            molecularweight += self.aaComp[aa]*(self.aa2mw[aa]-mwH20)




       #print (numberOfWater, sumOfMW)
        return molecularweight




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
