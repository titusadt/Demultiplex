#!/usr/bin/env python

import argparse
import bioinfo
import gzip
import itertools
import numpy as np
import matplotlib.pyplot as plt


def get_args():
    parser = argparse.ArgumentParser(description="A program for contig length distribution")
    parser.add_argument("-f1", "--file1", help="to specify the filename")#, required=True)
    parser.add_argument("-f2", "--file2", help="to specify the filename")#, required=True)
    parser.add_argument("-f3", "--file3", help="to specify the filename")#, required=True)
    parser.add_argument("-f4", "--file4", help="to specify the filename")#, required=True)
    parser.add_argument("-a", "--analysis_file", help="to specify the name of the analysis output file")#, required=True)
    parser.add_argument("-l", "--seq_length", help="to specify the legnth of the sequences", type=int)
    parser.add_argument("-p", "--plotname", help="to specify the name of the bar plot")#, required=True)
    return parser.parse_args()

args=get_args()
seq_length:int =args.seq_length
file1:str =args.file1
file2:str =args.file2
file3:str =args.file3
file4:str =args.file4
analysis:str=args.analysis_file
plotname:str=args.plotname


indexes='/projects/bgmp/shared/2017_sequencing/indexes.txt'



#creating counters to get the kown, unknow and blah blah blah countera
matched:int=0
hopped:int=0
unknown:int=0
total:int=0



def read_files(file)->tuple:
    '''This function takes in the filename and reads each line in the record and returns a tuple'''
    header = file.readline().strip()
    sequence = file.readline().strip()
    plus =  file.readline().strip()
    qual_score = file.readline().strip()
    return header, sequence, plus,qual_score

def rev_comp(index2):
    '''this function takes the second index (sequence in R3) and returns a reverse complement'''
    index2_reverse = index2[::-1]
    complement = {'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N'}
    index_rev_comp=[]
    for i in index2_reverse:
        index_rev_comp.append(complement[i])
    return ''.join(index_rev_comp)

def append_header(r1_header, r4_header, index1, rev_comp):
    '''this function takes the headers from r1 and r2 and append the index pairs to it and return the new header which will be written to the file'''
    index_pair= index1+'-'+rev_comp
    new_r1_header= r1_header+' '+index_pair
    new_r4_header = r4_header+' '+index_pair
    return new_r1_header, new_r4_header


#list to hold average quality score
qual_list:list=[]
def get_cutoff(qual_score):
    '''This function takes in the indexes and return True if any base is less than the cutoff and False if it is greater that the cutoff'''
    for seq in range(0,len(qual_score)):
        score=bioinfo.convert_phred(qual_score[seq])
    #check if the individual base pairs is less than the cutoff
        if score<26: 
            return True
       

#made the lookup index to create the output files
index_set=set()
with open(indexes, 'r') as index_file:
    index_file.readline()
    for line in index_file:
        line =line.split('\t')[4].strip('\n')
        index_set.add(line)


#matched_filename{key: indexs/barcodes (str), values:(filename_R1:str, filename_R2:str)}
matched_filename={}
def open_files():
    '''this function loops through the index sets and creates output files to write to'''
    for index in index_set:
        #create filenames to write to
        filename_R1 = index+'_R1.fastq'
        filename_R4 = index+'_R2.fastq'
        

        #open the filenames
        filename_R1 = open(filename_R1, 'w')
        filename_R4 = open(filename_R4, 'w')

        # print(filename_R1) 
        # print(filename_R4)

        #add them to matched_filename{key: indexex, values:(filename_R1, filename_R2)}
        file_tuple = (filename_R1, filename_R4)
        matched_filename[index]=file_tuple
    return matched_filename



unmatched_r1 = open('unknown_r1.fastq', 'w') 
unmatched_r2 = open('unknown_r2.fastq', 'w')
hopped_r1 =open('hopped_r1.fastq', 'w')
hopped_r2=open('hopped_r2.fastq', 'w')

matched_dict=open_files()
#print(matched_dict)

#key: index pairs, values: frequency of index pairs
matched_pairs:dict={}
#key: index pairs, values: frequency of the hopped pairs
hopped_pairs:dict={}
with gzip.open(file1, 'rt') as f1, gzip.open(file2, 'rt') as f2, gzip.open(file3, 'rt') as f3,gzip.open(file4, 'rt') as f4:
    
    while True:
        r1_header, r1_seq, r1_plus, r1_qual_score = read_files(f1)
        i1_header, i1_seq, i1_plus, i1_qual_score = read_files(f2)
        i2_header, i2_seq, i2_plus, i2_qual_score = read_files(f3)
        r4_header, r4_seq, r4_plus, r4_qual_score = read_files(f4)
        
        if r1_header == "" or r1_seq == "" or r1_plus == "" or r1_qual_score == "":
            break
        total+=1
        #get the reverse complement of the 2nd index
        reverse_comp = rev_comp(i2_seq)
        new_r1_header, new_r4_header=append_header(r1_header, r4_header, i1_seq, reverse_comp)
        index_pairs = i1_seq+'-'+reverse_comp
        #get the means for the indexes
        index1_mean = get_cutoff(i1_qual_score)
        index2_mean = get_cutoff(i2_qual_score)
        #print(f'index 1 mean: {index1_mean}\nindex 2 mean: {index2_mean}')
        
        #checking N is in index the indexes or if it is less than the cut off
        if 'N' in i1_seq or 'N' in reverse_comp: #or index1_mean ==True  or index2_mean ==True: #SPLIT THIS IF STATEMENT
            unknown+=1
            #print(unknown)
            #write unmatched R1 records to a new file, and write unmatched R2 records to a new file
            unmatched_r1.write(f'{new_r1_header}\n{r1_seq}\n{r1_plus}\n{r1_qual_score}\n')
            unmatched_r2.write(f'{new_r4_header}\n{r4_seq}\n{r4_plus}\n{r4_qual_score}\n')
        elif i1_seq not in matched_dict or reverse_comp not in matched_dict:
            unknown+=1
            unmatched_r1.write(f'{new_r1_header}\n{r1_seq}\n{r1_plus}\n{r1_qual_score}\n')
            unmatched_r2.write(f'{new_r4_header}\n{r4_seq}\n{r4_plus}\n{r4_qual_score}\n')
        #Now write the index matched files to their resperctive files
        #dictionaries ususally just loop through keys when no specification is made.
        elif i1_seq == reverse_comp: #and i1_seq in matched_dict:
            matched+=1
            if matched_pairs.get(index_pairs):
                matched_pairs[index_pairs]+=1
            else:
                matched_pairs[index_pairs]=1
            #for index, filename_tuple in matched_dict.items():
            r1_fh, r4_fh = matched_dict[i1_seq]
            r1_fh.write(f'{new_r1_header}\n{r1_seq}\n{r1_plus}\n{r1_qual_score}\n')
            r4_fh.write(f'{new_r4_header}\n{r4_seq}\n{r4_plus}\n{r4_qual_score}\n')
        elif i1_seq != reverse_comp: #and (i1_seq in matched_dict or reverse_comp in matched_dict): #check if theyre in the matched_dictionary the rest of them will go to unkown
            hopped+=1
            if hopped_pairs.get(index_pairs):
                hopped_pairs[index_pairs]+=1
            else:
                hopped_pairs[index_pairs]=1
            hopped_r1.write(f'{new_r1_header}\n{r1_seq}\n{r1_plus}\n{r1_qual_score}\n')
            hopped_r2.write(f'{new_r4_header}\n{r4_seq}\n{r4_plus}\n{r4_qual_score}\n')

#close the files youre writing too
unmatched_r1.close()
unmatched_r2.close()
hopped_r1.close()
hopped_r2.close()

#close the matched files
for r1_fh, r2_fh in matched_dict.values():
    r1_fh.close()
    r2_fh.close()


#plotting the matched pairs
plt.bar(matched_pairs.keys(), matched_pairs.values())
plt.title('Plot of Matcned pairs')
plt.xlabel('Matched Pairs')
plt.ylabel('Frequency')
plt.savefig(plotname)

#ANALYSIS AREA
with open(analysis, 'w') as output:
    #getting the total number of the matched unknown and hopped
    output.write(f'Number of Unknown: {unknown}\n')
    output.write(f'NUmber of Matched: {matched}\n')
    output.write(f'Number of hopped: {hopped}\n')
    output.write(f'Total amount of lines in the file: {total}\n')

    #Calculate the percentages 
    percent_unknown=(unknown/total)*100
    percent_matched=(matched/total)*100
    percent_hopped=(hopped/total)*100
    output.write(f'Percentage of unknown reads: {percent_unknown}\n')
    output.write(f'Percentage of matched reads: {percent_matched}\n')
    output.write(f'Percentage of hopped reads: {percent_hopped}\n')


   
    
#write out matched percentages to its own file
with open('matched_instances.txt', 'w') as matched_output:
    matched_output.write(f'Matched Pairs\tTotal\tTotal Percentage\tMatched Percentage\n')
    for key, value in matched_pairs.items():
        matched_total = (value/total)*100
        matched_percent = (value/matched)*100
        matched_output.write(f'{key}\t{value}\t{matched_total}\t{matched_percent}\n')

with open('hopped_instances.txt', 'w') as hopped_output:
    hopped_output.write(f'Hopped Pairs\tTotal\tTotal Percentage\tMatched Percentage\n')
    for key, value in hopped_pairs.items():
        hopped_total = (value/total)*100
        hopped_percent = (value/hopped)*100
        hopped_output.write(f'{key}\t{value}\t{hopped_total}\n{hopped_percent}\n')
