#!/usr/bin/env python

import bioinfo
import argparse
import gzip

def get_args():
    parser = argparse.ArgumentParser(description="This will demultiplex reads given 4 index files, a file of known index names and sequences /"
    "and a quality score cutoff.")
    parser.add_argument("-r1", "--read1", help="input read1 file name", type = str, required = True)
    parser.add_argument("-r2", "--read2", help="input read2 file name", type = str, required = True)
    parser.add_argument("-I1", "--index1", help="input index file name for read 1", type = str, required = True)
    parser.add_argument("-I2", "--index2", help="input index file name for read 2", type = str, required = True)
    parser.add_argument("-EI", "--expected_index_file", help="input the file that contains the list of expected index sequences and their names", type = str, required = True)
    parser.add_argument("-qc", "--quality_cutoff", help="input a number for the quality cutoff/threshold", type = int, required = True)
    return parser.parse_args()

#variables that hold the 4 read files, known index file, and the quality score cutoff
args = get_args()
Read1: str = args.read1
Read2: str = args.read2
Index1: str = args.index1
Index2: str = args.index2
known_index_file: str = args.expected_index_file
quality_cutoff = args.quality_cutoff

#initialize counters for the 3 conditions
total_matched = 0
total_hopped = 0
total_unknown = 0

#initialize dictionaries for each condition to hold the indexes as keys and their frequencies as values
hopped_dict: dict = {}
matched_dict: dict = {}

#make the dictionary for the known indices sequences (keys) and their names (values)
known_index_dict: dict = {}
#make the matched file names
matched_fn_dict: dict = {}

with open(known_index_file, "r") as indexes:
    for num, line in enumerate(indexes):
        if num != 0:
            line_list = line.split() #make a list of the items across each line
            index_name = line_list[3] #grab the index name
            index_seq = line_list[4] #grab the index sequence
            R1_file_name = f'{index_name}_R1.fq'
            R2_file_name = f'{index_name}_R2.fq'
            #store the index sequence as the key and the values as the index name
            known_index_dict[index_seq] = index_name
            #store the index name as the key and the values as a list of matched file names
            matched_fn_dict[f'{index_name}_R1'] = open(R1_file_name, "w")
            matched_fn_dict[f'{index_name}_R2'] = open(R2_file_name, "w")
            #initiate all indexes to value 0 in the matched dictionary that will count how many times that index was matching across all records
            matched_dict[index_name] = 0

#open all 52 files (matched, hopped, and unknown)
hopped_R1 = open("hopped_R1.fq", "w")
hopped_R2 = open("hopped_R2.fq", "w")
unknown_R1 = open("unknown_R1.fq", "w")
unknown_R2 = open("unknown_R2.fq", "w")

matched_fn_keys = list(matched_fn_dict.keys()) #put matched file name KEYS in a list

with gzip.open(Read1, 'rt') as R1, gzip.open(Read2, 'rt') as R2, gzip.open(Index1, 'rt') as I1, gzip.open(Index2, 'rt') as I2:
    
    while True:
        #go through each file and store the next record as it's parts (header, sequence, plus line, and qscore line)
        R1_header, R1_seq, R1_plus, R1_q_line = bioinfo.read_record(R1)
        R2_header, R2_seq, R2_plus, R2_q_line = bioinfo.read_record(R2)
        I1_header, I1_seq, I1_plus, I1_q_line = bioinfo.read_record(I1)
        I2_header, I2_seq, I2_plus, I2_q_line = bioinfo.read_record(I2)
        
        #break the while loop at the end of the file 
        if R1_header == "":
            break
        
        #update index 2 to the proper/expected index via reverse complement
        I2_seq = bioinfo.reverse_complement(I2_seq)

        #Prepare the header lines by appending them with the index pairs
        R1_header = f'{R1_header} {I1_seq}-{I2_seq}'
        R2_header = f'{R2_header} {I1_seq}-{I2_seq}'

        #Check for the 3 conditions (matched, hopped, and unknown)
        #Check if unknown (indexes are not in database or don't pass quality check)
        if I1_seq not in known_index_dict or I2_seq not in known_index_dict:
            total_unknown +=1 #increment total unkown counter
            #write out to the unknown R1 and R2 files 
            unknown_R1.write(f'{R1_header}\n{R1_seq}\n{R1_plus}\n{R1_q_line}\n')
            unknown_R2.write(f'{R2_header}\n{R2_seq}\n{R2_plus}\n{R2_q_line}\n')
        else: #everything that is has indexes in the dictionary
            #check if index qscore is good (more than 1 qscores for each position are less than cutoff for each index)
            I1_check: bool = bioinfo.qscore_check(I1_q_line, quality_cutoff, 1)
            I2_check: bool = True #initiate to qscore is good (True)
            if I1_check == True: #if index 1 passes quality, check index 2
                I2_check: bool = bioinfo.qscore_check(I2_q_line, quality_cutoff, 1)
                if I1_check == False or I2_check == False:
                    total_unknown +=1 #increment total unkown counter
                    #write out to the unknown R1 and R2 files 
                    unknown_R1.write(f'{R1_header}\n{R1_seq}\n{R1_plus}\n{R1_q_line}\n')
                    unknown_R2.write(f'{R2_header}\n{R2_seq}\n{R2_plus}\n{R2_q_line}\n')
                else: #everything that has indexes are in dictionary AND Indexes are good quality
                    #check if read qscore is good (more than 10% of bases have qscores less than the cutoff)
                    R1_check: bool = bioinfo.qscore_check(R1_q_line, quality_cutoff, int(len(R1_q_line)*0.1))
                    R2_check: bool = True
                    if R1_check == True: #if read 1 is good, check if read 2 is good
                        R2_check: bool = bioinfo.qscore_check(R2_q_line, quality_cutoff, int(len(R2_q_line)*0.1))
                        if R1_check == False or R2_check == False:
                            total_unknown +=1 #increment total unkown counter
                            #write out to the unknown R1 and R2 files 
                            unknown_R1.write(f'{R1_header}\n{R1_seq}\n{R1_plus}\n{R1_q_line}\n')
                            unknown_R2.write(f'{R2_header}\n{R2_seq}\n{R2_plus}\n{R2_q_line}\n')
                        else: #All indexes are in the dictionary AND all indexes and reads have good quality
                            if I1_seq == I2_seq: #check if matched
                                total_matched +=1 #increment total matched counter
                                index_name = known_index_dict[I1_seq] #grab the index name using the sequence
                                #increment match counter for this index
                                matched_dict[index_name] += 1
                                #write out to the matched R1 and R2 files
                                matched_fn_dict[f'{index_name}_R1'].write(f'{R1_header}\n{R1_seq}\n{R1_plus}\n{R1_q_line}\n')
                                matched_fn_dict[f'{index_name}_R2'].write(f'{R2_header}\n{R2_seq}\n{R2_plus}\n{R2_q_line}\n')
                            else: #must be "hopped"
                                total_hopped += 1 #increment the total hopped counter
                                index1_name = known_index_dict[I1_seq] #grab the index 1 name using the sequence
                                index2_name = known_index_dict[I2_seq] #grab the index 2 name using the sequence
                                #increment the counter for this index pair
                                if f'{index1_name}-{index2_name}' in hopped_dict: #is it already in the dictionary? -Yes
                                    hopped_dict[f'{index1_name}-{index2_name}'] +=1
                                if f'{index1_name}-{index2_name}' not in hopped_dict:  #it's not in the dictionary yet.
                                    hopped_dict[f'{index1_name}-{index2_name}'] = 1
                                #write out to the hopped R1 and R2 files
                                hopped_R1.write(f'{R1_header}\n{R1_seq}\n{R1_plus}\n{R1_q_line}\n')
                                hopped_R2.write(f'{R2_header}\n{R2_seq}\n{R2_plus}\n{R2_q_line}\n')

#Print Summary
#Matched Reads
print(f'Total Matched Reads: {total_matched}')
print("Matched Index Name   Number of Reads")
sorted_indexes = sorted(list(matched_dict.keys()))
for index in sorted_indexes:
    print(f'{index}\t{matched_dict[index]}')
#Hopped Reads
print(f'Total Hopped Reads: {total_hopped}')
print("Hopped Index Pairs   Number of Reads")
sorted_pairs = sorted(list(hopped_dict.keys()))
for index_pair in hopped_dict:
    print(f'{index_pair}\t{hopped_dict[index_pair]}')
#Unknown Reads
print(f'Total Unknown Reads: {total_unknown}')

#close all files
for pos, key in enumerate(matched_fn_keys):
    matched_fn_dict[key].close()

hopped_R1.close()
hopped_R2.close()
unknown_R1.close()
unknown_R2.close()







