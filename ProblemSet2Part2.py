
'''
program to “parse” sequence name information from a single line of a FASTQ formatted file from sequencer

Example:
converts this : @EAS139:136:FC706VJ:2:2104:15343:197393
to: 
Instrument = EAS139
Run ID = 136
Flow Cell ID = FC706VJ
Flow Cell Lane = 2
Tile Number = 2104
X-coord = 15343
Y-coord = 197393
'''


class FastqString(str):
    ''' Parses fast1 files'''

    def parse(self):
        ''' parses the input data'''
        #gets input
        Input = input('please print input here: ')
        #gets rid of the @ symbol
        stringWithoutAt = Input[1:]
        #splices string by :
        stringSplit = stringWithoutAt.split(":")
        return stringSplit



def main():
    ''' Function docstring goes here.'''
    #runs the input through the class objects
    fastq = FastqString(input)
    #runs the stringSplit throught class and object
    stringSplit = fastq.parse()
    #theses print lines print the coresponding splits values
    print("Instrument =", stringSplit[0])
    print("Run ID =", stringSplit[1])
    print("Flow Cell ID =", stringSplit[2])
    print("Flow Cell Lane =", stringSplit[3])
    print("Tile Number =", stringSplit[4])
    print("X-coord =", stringSplit[5])
    print("Y-coord =", stringSplit[6])


main()
