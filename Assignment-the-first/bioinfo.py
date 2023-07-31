#!/usr/bin/env python

# Author: <YOU> <optional@email.address>

# Check out some Python module resources:
#   - https://docs.python.org/3/tutorial/modules.html
#   - https://python101.pythonlibrary.org/chapter36_creating_modules_and_packages.html
#   - and many more: https://www.google.com/search?q=how+to+write+a+python+module

'''This module is a collection of useful bioinformatics functions
written during the Bioinformatics and Genomics Program coursework.
You should update this docstring to reflect what you would like it to say'''

__version__ = "0.5"         # Read way more about versioning here:
                            # https://en.wikipedia.org/wiki/Software_versioning

DNA_bases = set('ATGCatcg') 
RNA_bases = set('AUGCaucg') 

def convert_phred(letter: str) -> int:
    '''Converts a single character into a phred score'''
    #pass
    return ord(letter)-33

def qual_score(phred_score: str) -> float:
    '''Write your own doc string'''
    sum_score:int =0
    for single_letter in phred_score:
        sum_score+=convert_phred(single_letter)
    average=sum_score/len(phred_score)
    return average

def validate_base_seq(seq,RNAflag=False):
    '''This function takes a string. Returns True if string is composed
    of only As, Ts (or Us if RNAflag), Gs, Cs. False otherwise. Case insensitive.'''
    return set(seq)<=(RNA_bases if RNAflag else DNA_bases)


def gc_content(DNA):
    '''Returns GC content of a DNA or RNA sequence as a decimal between 0 and 1.'''
    assert validate_base_seq(DNA), "String contains invalid characters - are you sure you used a DNA sequence?"
    DNA = DNA.upper()
    return (DNA.count("G")+DNA.count("C"))/len(DNA)

# def oneline_fasta(filename, single_line):
#     '''This function takes in a fast file and writes out a new fast file with only single line sequence lines'''
#     with open(filename, 'r') as fa_file, open(single_line, 'w') as outfile:
#     sequence: str=''
#     first_line=True
#     for line in fa_file:
#         line =line.strip()
#         if first_line:
#             header = line
#             outfile.write(f'{header}\n')
#             first_line=False
#         elif not line.startswith('>'):
#             sequence+=line.strip('\n')
#         elif line.startswith('>'):
#             header2=line
#             outfile.write(f'{sequence}\n')
#             outfile.write(f'{header2}\n')
#             sequence=""
#     outfile.write(f'{sequence}\n')

def calc_median(sort_qscores): #-> tuple[list, int]:\n",
    """This function calculates the median of each line of quality scores"""
    #check if the line is even
    length = len(sort_qscores)
    if length%2==0:
        mid1 = sort_qscores[length//2]
        mid2 = sort_qscores[length//2-1]
        median = (mid1+mid2)/2
    else:
        median = sort_qscores[length//2]
    return median


if __name__ == "__main__":
    # write tests for functions above, Leslie has already populated some tests for convert_phred
    assert convert_phred("I") == 40, "wrong phred score for 'I'"
    assert convert_phred("C") == 34, "wrong phred score for 'C'"
    assert convert_phred("2") == 17, "wrong phred score for '2'"
    assert convert_phred("@") == 31, "wrong phred score for '@'"
    assert convert_phred("$") == 3, "wrong phred score for '$'"
    print("Your convert_phred function is working! Nice job")

    #tests for validate_base_seq
    print("testing RNA:")
    assert validate_base_seq("AAUAGAU", True), "RNA test failed"
    print("RNA test passed!")
    print("testing DNA:")
    assert validate_base_seq("AATAGAT"), "DNA test failed"
    print("DNA test passed!")
    print("testing non-nucleic acid:")
    assert validate_base_seq("R is the best!")==False, "R sux"
    print("non-nucleic test passsed!")

    #test for qual score
    assert qual_score("EEE") == 36
    assert qual_score("#I") == 21
    assert qual_score("EJ") == 38.5
    #assert qual_score(37.62105263157895) == 37.62105263157895, "wrong average phred score"
    print("You calcluated the correct average phred score")
