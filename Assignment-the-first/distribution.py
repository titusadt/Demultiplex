#!/usr/bin/env python
import bioinfo
import argparse
import matplotlib.pyplot as plt
import gzip



#GONNA ARGPARSE THE FILE INTO THIS
def get_args():
    parser = argparse.ArgumentParser(description="A program for contig length distribution")
    parser.add_argument("-f1", "--file1", help="to specify the filename")#, required=True)
    parser.add_argument("-l", "--seq_length", help="to specify the legnth of the sequences", type=int)
    parser.add_argument("-p", "--plotname", help="to specify the name of the bar plot")#, required=True)
    return parser.parse_args()

#takes on the values from get_args
args=get_args()
file1:str =args.file1
seq_length:int =args.seq_length
plotname:str = args.plotname

def populate_list(file1) -> tuple[list, int]:
    """This function populates an empty list with their phred scores"""
    empty_list: list = []
    empty_list= [0]*seq_length
    count: int =0
    with gzip.open(file1, 'rt') as file1:
        i: int =0
        for line in file1:
            i+=1
            count+=1
            #print(line[0])
            line=line.strip('\n')
            if i%4 ==0:
                for seq in range(0,len(line)):
                    empty_list[seq]+=bioinfo.convert_phred(line[seq])
        
    return empty_list, count

#running 
my_list, num_lines = populate_list(file1)

for base in range(len(my_list)):
    my_list[base]=my_list[base]/(num_lines/4) #(len(my_list))\n",



plt.bar(range(seq_length), my_list)
plt.xlabel('Mean Score')
plt.ylabel('Base Pair')
plt.title('Mean score at each base pair')
plt.savefig(plotname)
