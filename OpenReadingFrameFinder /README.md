This program analyizes open reading frame of DNA sequences based on Start Codon Sequence 'AUG' and a minimum gene length of 100bp 

This program also makes the assumption that if there is a stop codon and no start codon from the left then it will use the left end of the genome as the start point

Requires FASTA File 
A file has been provided in the folder called 'testfa.fa'
The program has this as the input, to change to standard in please change main

example output: 

(name of sequence)

(reading Frame) (location on sequence) (length of reading frame)

tass2 NODE_159_length_75728_cov_97.549133

+1 57169..61908 4740

+1 57184..61908 4725

+1 57190..61908 4719

+1 57235..61908 4674

+1 57319..61908 4590

