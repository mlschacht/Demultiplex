#!/usr/bin/env python

import bioinfo
import argparse

#argparse CHANGE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# def get_args():
#     parser = argparse.ArgumentParser(description="This will find the longest protein for each gene /"
#     "and output it to its own fastq file.")
#     parser.add_argument("-f", "--filename", help="input file name", type = str, required = True)
#     parser.add_argument("-o", "--outputfilename", help="what do you want the output file to be named?", type = str, required = True)
#     return parser.parse_args()

# args = get_args()
# fn: str = args.filename #file name
# outfile: str = args.outputfilename #desired output file name

#variables that hold the 4 read files
#update after adding argparse
Read1 = "/projects/bgmp/mlscha/bioinfo/Bi622/Demultiplex/TEST-input_FASTQ/R1_input_test.fq"
Read2 = "/projects/bgmp/mlscha/bioinfo/Bi622/Demultiplex/TEST-input_FASTQ/R4_input_test.fq"
Index1 = "/projects/bgmp/mlscha/bioinfo/Bi622/Demultiplex/TEST-input_FASTQ/R2_input_test.fq"
Index2 = "/projects/bgmp/mlscha/bioinfo/Bi622/Demultiplex/TEST-input_FASTQ/R3_input_test.fq"
known_index_file = "/projects/bgmp/shared/2017_sequencing/indexes.txt"
quality_cutoff = 29

#initialize counters for the 3 conditions
total_matched = 0
total_hopped = 0
total_unknown = 0

#initialize dictionaries for each condition to hold the indexes as keys and their frequencies as values
#unknown_dict: dict = {}
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

# for pos, key in enumerate(matched_fn_keys): #MIGHT NOT NEED THIS ANYMORE
#     file_name = matched_fn_dict[key]
#     matched_fn_dict[key] = open(file_name, "w") #open each of the matched files

with open(Read1, 'r') as R1, open(Read2, 'r') as R2, open(Index1, 'r') as I1, open(Index2, 'r') as I2:
    
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

        #check if index qscore is good (more than 1 qscores for each position are less than cutoff for each index)
        I1_check: bool = bioinfo.qscore_check(I1_q_line, quality_cutoff, 1)
        I2_check: bool = True #initiate to qscore is good (True)
        if I1_check == True: #if index 1 passes quality, check index 2
            I2_check: bool = bioinfo.qscore_check(I2_q_line, quality_cutoff, 1)

        #check if read qscore is good (more than 10% of bases have qscores less than the cutoff)
        print(int(len(R1_q_line)*0.1))
        #R1_check: bool = bioinfo.qscore_check(R1_q_line, quality_cutoff, int(len(R1_q_line)*0.1)

        #Prepare the header lines by appending them with the index pairs
        R1_header = f'{R1_header} {I1_seq}-{I2_seq}'
        R2_header = f'{R2_header} {I1_seq}-{I2_seq}'

        #Check for the 3 conditions (matched, hopped, and unknown)
        #Check if unknown (indexes are not in database or don't pass quality check)
        if I1_seq not in known_index_dict or I2_seq not in known_index_dict or I1_check == False or I2_check == False:
            total_unknown +=1 #increment total unkown counter
            #write out to the unknown R1 and R2 files 
            if total_unknown > 1: #check if this is not the first record (add new line at the beginning)
                unknown_R1.write(f'\n{R1_header}\n{R1_seq}\n{R1_plus}\n{R1_q_line}')
                unknown_R2.write(f'\n{R2_header}\n{R2_seq}\n{R2_plus}\n{R2_q_line}')
            else: #if it is the first record, don't add a new line at the beginning
                unknown_R1.write(f'{R1_header}\n{R1_seq}\n{R1_plus}\n{R1_q_line}')
                unknown_R2.write(f'{R2_header}\n{R2_seq}\n{R2_plus}\n{R2_q_line}')
        elif I1_seq == I2_seq: #check if matched
            total_matched +=1 #increment total matched counter
            index_name = known_index_dict[I1_seq] #grab the index name using the sequence
            #increment match counter for this index
            matched_dict[index_name] += 1
            #write out to the matched R1 and R2 files
            if total_matched >1: #check if this is not the first record (add new line at the beginning)
                matched_fn_dict[f'{index_name}_R1'].write(f'\n{R1_header}\n{R1_seq}\n{R1_plus}\n{R1_q_line}')
                matched_fn_dict[f'{index_name}_R2'].write(f'\n{R2_header}\n{R2_seq}\n{R2_plus}\n{R2_q_line}')
            else: #if it is the first record, don't add a new line at the beginning
                matched_fn_dict[f'{index_name}_R1'].write(f'{R1_header}\n{R1_seq}\n{R1_plus}\n{R1_q_line}')
                matched_fn_dict[f'{index_name}_R2'].write(f'{R2_header}\n{R2_seq}\n{R2_plus}\n{R2_q_line}')
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
            if total_hopped > 1: #check if this is not the first record (add new line at the beginning)
                hopped_R1.write(f'\n{R1_header}\n{R1_seq}\n{R1_plus}\n{R1_q_line}')
                hopped_R2.write(f'\n{R2_header}\n{R2_seq}\n{R2_plus}\n{R2_q_line}')
            else: #if it is the first record, don't add a new line at the beginning
                hopped_R1.write(f'{R1_header}\n{R1_seq}\n{R1_plus}\n{R1_q_line}')
                hopped_R2.write(f'{R2_header}\n{R2_seq}\n{R2_plus}\n{R2_q_line}')

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







