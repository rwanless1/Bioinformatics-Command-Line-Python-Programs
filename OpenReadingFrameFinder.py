#!/usr/bin/env python3
# Name: Ryan Wanless (rwanless)
# Group Members: Boxuan Ma(bma79), Jesse Smith(jefsmith), Shreya Sinha(ssinha10)
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

from FastaReader import FastAreader

########################################################################
# CommandLine
########################################################################
class CommandLine():
    '''
    Handle the command line, usage and help requests.

    CommandLine uses argparse, now standard in 2.7 and beyond.
    it implements a standard command line argument parser with various argument options,
    a standard usage and help.

    attributes:
    all arguments received from the commandline using .add_argument will be
    avalable within the .args attribute of object instantiated from CommandLine.
    For example, if myCommandLine is an object of the class, and requiredbool was
    set as an option using add_argument, then myCommandLine.args.requiredbool will
    name that option.

    '''

    def __init__(self, inOpts=None):
        '''
        Implement a parser to interpret the command line argv string using argparse.
        '''

        import argparse
        self.parser = argparse.ArgumentParser(
            description='Program prolog - a brief description of what this thing does',
            epilog='Program epilog - some other stuff you feel compelled to say',
            add_help=True,  # default is True
            prefix_chars='-',
            usage='%(prog)s [options] -option1[default] <input >output'
            )
        self.parser.add_argument('-lG', '--longestGene', action='store', nargs='?', const=True, default=False,
                                 help='longest Gene in an ORF')
        self.parser.add_argument('-mG', '--minGene', type=int, choices=(100, 200, 300, 500, 1000), default=100,
                                 action='store', help='minimum Gene length')
        self.parser.add_argument('-s', '--start', action='append', default=['ATG'], nargs='?',
                                 help='start Codon')  # allows multiple list options
        self.parser.add_argument('-t', '--stop', action='append', default=['TAG', 'TGA', 'TAA'], nargs='?',
                                 help='stop Codon')  # allows multiple list options
        self.parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.1')
        if inOpts is None:
            self.args = self.parser.parse_args()
        else:
            self.args = self.parser.parse_args(inOpts)


'''Used for finding the OFS '''
class Codonsearch():
    '''constructor of the sequence with the stop codons and the sequence'''
    def __init__(self, seq):

        # creates the list of dictionaries which are of the individual codon frames
        self.seq = seq
        self.StopCodons = ['TAG', 'TGA', 'TAA']





    '''This will find 6 codon frames'''
    def CodonFind(self, minimum = 100, longestGene = False ):
        self.seqlen = len(self.seq)
        self.StopCodons = ['TAG', 'TGA', 'TAA']
        self.geneCanadates = []
        self.startPosistionList = []
        self.startcodons = set(('ATG', 'GTG', 'TTG'))


        #used for loop with the range of the sequence and then iterates for every three codons which is a frame
        #start position with [0,] handless the edge case
        for frame in range(3):
            self.startPosistionList = [0,]
            for i in range(frame, len(self.seq), 3):
                codon = self.seq[i: i+3]

                if (codon in self.startcodons and longestGene == False) or (longestGene == True and len(self.startPosistionList) < 1):
                    self.startPosistionList.append(i)
                if codon in self.StopCodons:
                    for start in self.startPosistionList:
                        length = (i + 3) - (start + 3)
                        if length > minimum:
                            #appends the list to show frame then start position then end, then the length of the orf
                            self.geneCanadates.append([frame + 1, (start +4), i+3, ((i+3) - (start + 3))])
                    self.startPosistionList = []
        for start in self.startPosistionList:
            length = (i + 3) - (start + 3)
            if length > minimum:
               self.geneCanadates.append([frame + 1, (start + 4) , i+3, ((i+3) - (start + 3))])

        #reverse compliment of the sequence to handle frames 4 to 6

        self.reverseCompliment = self.seq.replace('A', 't').replace('T', 'a').replace('G', 'c').replace('C','g')
        self.reverseCompliment = self.reverseCompliment.upper()
        self.reverseCompliment = self.reverseCompliment[::-1]

        # does the same as the first for loop except does the reverse compliment
        for frame in range(3):
            self.startPosistionList = [0,]
            for i in range(frame, len(self.seq), 3):
                codon = self.reverseCompliment[i: i+3]

                if (codon in self.startcodons and longestGene == False) or (longestGene == True and len(self.startPosistionList) < 1):
                    self.startPosistionList.append(i)
                if codon in self.StopCodons:
                    for start in self.startPosistionList:
                        length = (i+3) - (start + 3)
                        if length > minimum:
                            self.geneCanadates.append([-(frame + 1), start, i + 3, ((i + 3) - (start + 3))])
                    self.startPosistionList = []
        for start in self.startPosistionList:
            length = (i + 3) - (start + 3)
            if length > minimum:
                self.geneCanadates.append([ -(frame +1), start, i+3, ((i+3) - (start + 3))])

    '''returns the orfs that was found in the sequence'''
    def geneCanad(self):
        return self.geneCanadates


def main(inFile=None, options=None):
    '''
    Find some genes.
    '''

    thisCommandLine = CommandLine(options)
    reader = FastAreader(inFile)
    for head, seq in reader.readFasta():
        orf = Codonsearch (seq)
        orf.CodonFind(thisCommandLine.args.minGene, thisCommandLine.args.longestGene)
        orf.geneCanad()
        print(head)
        #this line sorts the geneCanates list to the input that was required
        sortedGeneCandates = orf.geneCanad()
        sortedGeneCandates.sort(key=lambda x:(-x[3], x[1]))
        #sorts the gene canadits list based on the predetermined output
        for canidate in sortedGeneCandates:
            print('{:+d} {:>5d}..{:>5d} {:>d}'.format(canidate[0], canidate[1], canidate[2], canidate[3]))








if __name__ == "__main__":
    main('testfa.fa')  # delete this stuff if running from commandline


