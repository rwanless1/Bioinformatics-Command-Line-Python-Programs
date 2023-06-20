

'''
“clean up” a sequence of DNA by removing ambiguous bases (denoted by “N”) output from a sequencer.
'''
class DNAstring(str):
    def length(self):
        return (length(self))

    def purify(self):
        ''' Return an upcased version of the string, collapsing a single run of Ns.'''
        pureDNA = self.upper()
        # gets the left string by getting the correspoding values to the left of N
        left = pureDNA[:pureDNA.find('N')]
        # counts the number of Ns
        Nblock = pureDNA.count('N')
        # gets the right by finding length of the left plus the Ns and finging everything after that
        right = pureDNA[len(left) + Nblock:]
        # prints the left string, the right string, and the N numbers togehte randn formates them
        print(left + "{", Nblock, '}' + right)

        return pureDNA


def main():
    # Get user DNA data and clean it up.
    while (True):
        data = input('DNA data?')
        thisDNA = DNAstring(data)
        pureData = thisDNA.purify()


main()
