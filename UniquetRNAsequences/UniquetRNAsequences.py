

import sys


class FastAreader:

    def __init__(self, fname=''):
        '''contructor: saves attribute fname '''
        self.fname = fname

    def doOpen(self):
        if self.fname == '':
            return sys.stdin
        else:
            return open(self.fname)

    def readFasta(self):

        header = ''
        sequence = ''

        with self.doOpen() as fileH:

            header = ''
            sequence = ''

            # skip to first fasta header
            line = fileH.readline()
            while not line.startswith('>'):
                if not line:  # we are at EOF
                    return header, sequence
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
        pass

class tRNA:
    tRNAList= []
    '''Creates the lists and sequence input for the later objects'''
    def __init__(self, header, sequence):
        self.headers = header # creates a variable with the header
        self.sequences = sequence # creates sequennce variable
        self.powerSet = self.getPowerSet() #creates variable that calles getPowerSet method
        self.unique = set()
        self.essentialSet = set()
        tRNA.tRNAList.append(self) # creates a list of all tRNAs

    '''creates the powerset of the input sequence'''
    def getPowerSet(self):
        powerSet = set()
        for i in range(len(self.sequences)):#loops though each nucleotide
            for j in range(i+1,len(self.sequences)): # goes through the next nucleotide in the sequence
                substring = self.sequences[i:j] # this creates all the possible substrings of the sequences
                powerSet.add(substring) # adds all the possible substrings to the powerset variable
        return powerSet


    '''creates the set of unique charaters'''
    def getUnique(self):
        uniqueSet = set()
        for thisTRNA in tRNA.tRNAList:
            if thisTRNA is not self: # checks to see if this subset is already present in self
                duplicateSet = uniqueSet.union(thisTRNA.powerSet) #adds the duplicates to variable to duplicate set
        self.unique = self.powerSet - duplicateSet #subtracts the duplicates from the powerset to ge the uniques
        return self.unique


    '''Creates the list of essential tRNAs'''
    def getessential(self):
        nonEssential = set()
        for unique in self.unique:#loops through and eliminates ends of the sequence
            left = unique[:-1]
            right = unique[1:]

            if left in self.unique: #checks to see if the cut sequences are in unique and adds them to nonessential if they are
                nonEssential.add(unique)
            elif right in self.unique:
                nonEssential.add(unique)
        self.essentialSet = self.unique - nonEssential #creates the essential set by subtracting unique minus essential

        return self.essentialSet






    '''prints and formats the output'''
    def outPut(self):



        for essentialSetIndex in range(len(self.essentialList)): #loops through the essential list
            individualSequenceIndex = []
            currentSeq = self.sequences[essentialSetIndex]
            for essentialString in self.essentialList[essentialSetIndex]:
                currentIndex = currentSeq.find(essentialString)
                individualSequenceIndex.append([essentialString, currentIndex]) #finds the length of the essential item from the sequence

            individualSequenceIndex.sort(key=lambda item:item[1]) #prints out the output
            for lits in individualSequenceIndex:
                print('.' * lits[1] + lits[0])


def main(inFile = 'bos-tRNA.fa'):
    myFasta = FastAreader(inFile)


    for head, seq in myFasta.readFasta():
        cleanSeq = seq.replace('-', '').replace('_', '',).replace(',','')
        tRNA(head, cleanSeq)

    for thisTRNa in tRNA.tRNAList:
        thisTRNa.getUnique()
        thisTRNa.getessential()
    tRNA.tRNAList.sort(key= lambda x: x.headers)
    for thisTRNa in tRNA.tRNAList:
        print(thisTRNa.headers)
        print(thisTRNa.sequences)
        output = []
        for essential in thisTRNa.essentialSet:
            output.append('.' * thisTRNa.sequences.find(essential) + essential)
        output.sort(key=len)
        for item in output:
            print(item)





if __name__ == "__main__":
  main('bos-tRNA.fa')
