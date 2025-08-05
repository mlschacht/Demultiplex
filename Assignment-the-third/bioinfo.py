#!/usr/bin/env python
# Author: Makayla <mlscha@uoregon.edu>

# Check out some Python module resources:
#   - https://docs.python.org/3/tutorial/modules.html
#   - https://python101.pythonlibrary.org/chapter36_creating_modules_and_packages.html
#   - and many more: https://www.google.com/search?q=how+to+write+a+python+module

'''Helpful Bioinformatics Functions :)'''

__version__ = "5"         # Read way more about versioning here:
                            # https://en.wikipedia.org/wiki/Software_versioning
import gzip

DNA_bases = "ACGTNacgtn"
RNA_bases = "ACGUNacgun"
reverse_comp_dict: dict = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N':'N'}

def convert_phred(letter: str) -> int: #from PS3
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

def validate_base_seq(seq, RNAflag=False): #from Python 5 Notes
    '''This function takes a string. Returns True if string is composed
    of only As, Ts (or Us if RNAflag), Gs, Cs. False otherwise. Case insensitive.'''
    valid_bases = RNA_bases if RNAflag else DNA_bases
    return all([base in valid_bases for base in seq.upper()]) 

def gc_content(seq, RNAflag=False): #from Python notes 2
    '''Returns GC content of a DNA or RNA sequence as a decimal between 0 and 1, \
    unless a valid DNA/RNA sequence was not given'''
    assert validate_base_seq(seq, RNAflag) == True, \
     "This is not a valid sequence, or this is RNA and the RNAflag needs to be set to True"
    seq = seq.upper()
    gc = seq.count("G") + seq.count("C")
    num = gc / len(seq)
    return num

def calc_median(lst): #from PS4
    '''Given a list, returns the median value of the list'''
    lst.sort()
    half_way = len(lst) // 2
    if len(lst) % 2 == 0:
        median = (lst[half_way] + lst[half_way-1]) /2
    else:
        median = lst[half_way]
    return median

def oneline_fasta(multi_fa: str, single_fa: str)-> None:
    '''This function takes a multilined fasta and the name you want for your new file'''
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

def reverse_complement(seq: str, RNAflag=False) -> str:
    assert validate_base_seq(seq, RNAflag) == True, \
     "This is not a valid sequence, or this is RNA and the RNAflag needs to be set to True"
    rev_comp = ""
    for pos, base in enumerate(reversed(seq)):
        rev_comp = f'{rev_comp}{reverse_comp_dict[base]}'
    return rev_comp

def read_record(file) -> list:
    """Reads through each record of the file and returns that record."""
    record = [file.readline().strip(), file.readline().strip(), file.readline().strip(), file.readline().strip()]
    return record

def qscore_check(index: str, cutoff: int, max: int) -> bool:
    """checks if the index has a good qscore (True) or not a good qscore (False) more than 1 quality score below the cutoff"""
    counter = 0
    check = True #initiate to yes, this is a good qscore
    for i in index:
        q = convert_phred(i)
        if q < cutoff:
            counter += 1
        if counter > max:
            check = False #more than 1 is bad, this is not a good qscore
            break
    return check

if __name__ == "__main__":
    # write tests for functions above, Leslie has already populated some tests for convert_phred
    # These tests are run when you execute this file directly (instead of importing it)
    assert convert_phred("I") == 40, "wrong phred score for 'I'"
    assert convert_phred("C") == 34, "wrong phred score for 'C'"
    assert convert_phred("2") == 17, "wrong phred score for '2'"
    assert convert_phred("@") == 31, "wrong phred score for '@'"
    assert convert_phred("$") == 3, "wrong phred score for '$'"
    print("Your convert_phred function is working! Nice job")
    assert gc_content("GCGCGC") == 1, "wrong gc content"
    assert gc_content("AATTATA") == 0, "wrong gc content"
    assert gc_content("GCATCGAT") == 0.5, "wrong gc content"
    assert gc_content("acgtacgt") == 0.5, "can it handle lowercase strings?"
    print("correctly calculated GC content")
    assert qual_score("ABCDE") == 34, "wrong average phred score"
    assert qual_score("EEE") == 36, "wrong average phred score"
    assert qual_score("#I") == 21, "wrong average phred score"
    assert qual_score("EJ") == 38.5, "wrong average phred score"
    assert qual_score('FFHHHHHJJJJIJIJJJIJJJJJJIIIJJJEHJJJJJJJIJIDGEHIJJFIGGGHFGHGFFF@EEDE@C??DDDDDDD@CDDDDBBDDDBDBDD@') == 37.62105263157895, "wrong average phred score"
    print("Your qual_score function is working! Nice job")
    assert calc_median([1, 1, 1, 1]) == 1, "wrong median"
    assert calc_median([0, 1, 2, 3, 4, 5, 6, 7, 8]) == 4, "wrong median"
    assert calc_median([0, 8, 2, 4, 5]) == 4, "can it handle an unsorted list?"
    print("Your calc_median function is working! Nice job")
    assert validate_base_seq("AAAAATTTT", False), "Does it recognize A's, T's and DNA?"
    assert validate_base_seq("AAAAAUUUU", True), "Does it recognize A's, U's and RNA?"
    assert validate_base_seq("acgtactg", False), "Does it recognize lowercase?"
    assert validate_base_seq("acgtactgBBBB", False) == False, "Does it recognize non-DNA/RNA sequences?"
    print("Your validate_base_sequence function is working! Nice job")
    assert reverse_complement("ACT") == "AGT", "wrong reverse complement for 'ACT'"
    print("Your reverse complement function is working! Nice job")
    assert qscore_check("IIII", 30, 1) == True, "This is a good quality score..."
    assert qscore_check("###", 30, 1) == False, "This is a BAD quality score..."
    assert qscore_check("IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII###########", 29, 10) == False, "This is a BAD quality score..."
    print("Your qscore check function is working! Nice job")

    #gc_content("not DNA nor RNA string")  #This should fail
    #reverse_complement("not DNA nor RNA string")  #This should fail