#!/usr/bin/env python

# Author: <YOU> <optional@email.address>

# Check out some Python module resources:
#   - https://docs.python.org/3/tutorial/modules.html
#   - https://python101.pythonlibrary.org/chapter36_creating_modules_and_packages.html
#   - and many more: https://www.google.com/search?q=how+to+write+a+python+module

'''This module is a collection of useful bioinformatics functions
written during the Bioinformatics and Genomics Program coursework.
You should update this docstring to reflect what you would like it to say'''

__version__ = "5"         # Read way more about versioning here:
                            # https://en.wikipedia.org/wiki/Software_versioning

DNA_bases = "ACGTNacgtn"
RNA_bases = "ACGUNacgun"

import gzip

def convert_phred(letter: str) -> int:
    '''Converts a single character into a phred score'''
    return ord(letter) - 33

def qual_score(phred_score: str) -> float: #from PS3
    """takes the original, unmodified phred_score string as a parameter. 
    This function should calculate the average quality score of the whole phred string. 
    Be sure to write a doc string for this funcion."""
    #set total to 0
    total = 0
    # convert phred score
    for char in phred_score:
        q = convert_phred(char)
        #add it to the total
        total +=q
    #take the average
    avg = total / len(phred_score)
    #return the average
    return avg

def validate_base_sequence(seq, RNAflag=False): #from Python 5 Notes, this one needs checked???
    '''This function takes a string. Returns True if string is composed
    of only As, Ts (or Us if RNAflag), Gs, Cs. False otherwise. Case insensitive.'''
    valid_bases = RNA_bases if RNAflag else DNA_bases
    return all([base in valid_bases for base in seq.upper()])

def gc_content():
    '''Returns GC content of a DNA or RNA sequence as a decimal between 0 and 1.'''
    pass

def calc_median(lst): #from PS4, don't I need to have arguments??
    '''Given a sorted list, returns the median value of the list'''
    half_way = len(lst) // 2
    if len(lst) % 2 == 0:
        median = (lst[half_way] + lst[half_way-1]) /2
    else:
        median = lst[half_way]
    return median

def oneline_fasta(multi_fa: str, single_fa: str)-> None:
    '''This function takes a file and the name you want for your new file'''
    #The only function where checks are not needed/required.
    with open(multi_fa, 'r') as multi, open(single_fa, 'w') as single:
        for i, line in enumerate(multi):
            if line[0] == '>':
                if i!=0:
                    single.write('\n')
                single.write(line)
            else:
                single.write(line.strip('\n'))

        single.write("\n")
    return None

def init_list(lst: list, value: float=0.0) -> list:
    '''This function takes an empty list and will populate it with
    the value passed in "value". If no value is passed, initializes list
    with 101 values of 0.0.'''
    i=0
    while i < 101:
        lst.append(value)
        i+=1
    return lst

def populate_list(file: str) -> tuple[list, int]:
    """Loops through each line in the file to add up the phred scores at each position.
       Stores the total phred scores at each position in a list.
       Counts the total num of lines in the file.
       Returns the list of total phred scores and total line count."""
    qscore_list = []
    qscore_list = init_list(qscore_list, 0.0) #creates an empty list called qscore_list
    i = 0 #num lines
    with open(file, "r") as fq: #opens fastq file inputted
        for line in fq: #go through each line
            line = line.strip()
            if i%4==3: #check if this is a qscore line, not anything else
                for q, char in enumerate(line): #and go through each char in the line
                    phred = convert_phred(char) #convert the char to phred score
                    qscore_list[q] += phred #"add" the phred score to ongoing sum inside quality score list
            i+=1 #go to next line
            
    return qscore_list, i

def populate_list_gzip(file: str) -> tuple[list, int]:
    """Loops through each line in the file to add up the phred scores at each position.
       Stores the total phred scores at each position in a list.
       Counts the total num of lines in the file.
       Returns the list of total phred scores and total line count."""
    qscore_list = []
    qscore_list = init_list(qscore_list, 0.0) #creates an empty list called qscore_list
    i = 0 #num lines
    with gzip.open(file, "rt") as fq: #opens fastq file inputted
        for line in fq: #go through each line
            line = line.strip()
            if i%4==3: #check if this is a qscore line, not anything else
                for q, char in enumerate(line): #and go through each char in the line
                    phred = convert_phred(char) #convert the char to phred score
                    qscore_list[q] += phred #"add" the phred score to ongoing sum inside quality score list
            i+=1 #go to next line
            
    return qscore_list, i

if __name__ == "__main__":
    # write tests for functions above, Leslie has already populated some tests for convert_phred
    # These tests are run when you execute this file directly (instead of importing it)
    assert convert_phred("I") == 40, "wrong phred score for 'I'"
    assert convert_phred("C") == 34, "wrong phred score for 'C'"
    assert convert_phred("2") == 17, "wrong phred score for '2'"
    assert convert_phred("@") == 31, "wrong phred score for '@'"
    assert convert_phred("$") == 3, "wrong phred score for '$'"
    print("Your convert_phred function is working! Nice job")
