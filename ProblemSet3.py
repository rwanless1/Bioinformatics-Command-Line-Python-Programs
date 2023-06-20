"This class will parse through an amino acid sequencec and out put molecular weight, aacomposition, charge, PI"
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

    def __init__(self, protein):
        '''this is used for creating the dictionary which is used by the later objects to store the proper amino acids which 
        which are actually here.'''

        self.aaComp = {
            'A': 0, 'G': 0, 'M': 0, 'S': 0, 'C': 0,
            'H': 0, 'N': 0, 'T': 0, 'D': 0, 'I': 0,
            'P': 0, 'V': 0, 'E': 0, 'K': 0, 'Q': 0,
            'W': 0, 'F': 0, 'L': 0, 'R': 0, 'Y': 0
            }
            # count symbols in protein
            # ignoring any bad characters
        for aa in protein.upper():
            if aa in self.aaComp:  # count valid AA
                self.aaComp[aa] += 1
    'counts amino acid values'
    def aaCount(self):
        self.aaCountV = sum (self.aaComp.values()) #counts values from aaComp table and saves them to variable
        return self.aaCountV
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
        if self.aaCountV > 0:
            for key in self.aaComp.keys():
                percentage = (self.aaComp[key] / self.aaCount()*100) # then multiplies by 100 to get the percent and 0.2f for two decimal places
                print(key, '=','%.2f' % percentage,'%' )
        else:
            for key in self.aaComp.keys():
                percentage = (self.aaComp[key] / 1) # then multiplies by 100 to get the percent and 0.2f for two decimal places
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





        return molecularweight


# Please do not modify any of the following.  This will produce a standard output that can be parsed

import sys


def main():
    inString = input('protein sequence?')
    while inString:
        myParamMaker = ProteinParam(inString)
        myAAnumber = myParamMaker.aaCount()
        print("Number of Amino Acids: {aaNum}".format(aaNum=myAAnumber))
        print("Molecular Weight: {:.1f}".format(myParamMaker.molecularWeight()))
        print("molar Extinction coefficient: {:.2f}".format(myParamMaker.molarExtinction()))
        print("mass Extinction coefficient: {:.2f}".format(myParamMaker.massExtinction()))
        print("Theoretical pI: {ans}".format(ans=myParamMaker.pI()))
        print("Amino acid composition:")
        myParamMaker.aaComposition()

       # if myAAnumber == 0: myAAnumber = 1  # handles the case where no AA are present
        #print(myParamMaker.aaComposition())

        #for aa, n in sorted(myParamMaker.aaComposition(),
                              # key=lambda item: item[0]):
           # print("\t{} = {:.2%}".format(aa, n / myAAnumber))
        break
        self.protein = input('protein sequence?')


if __name__ == "__main__":
    main()
