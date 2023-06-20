
'''
program to “parse” sequence name information from a single line of a FASTQ formatted file
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
