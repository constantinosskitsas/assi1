########################################################
##Read in files
########################################################

def parse_fasta(path):
    fasta_list = {}
    fasta_file = open(path,"r")
    fasta = fasta_file.read()
    total = fasta.split(">")[1:]
    for i in total:
        name_seq = i.split()
        name = name_seq[0]
        seq = "".join(name_seq[1:])
        fasta_list[name] = seq
        
    fasta_file.close()
    return fasta_list


def parse_fastq(path):
    fastq_list = {}
    fastq_file = open(path,"r")
    fastq = fastq_file.read()
    reads=fastq.split("@")[1:]
    for read in reads:
        nameseq_qual = read.split("+")
        name_seq = nameseq_qual[0].split()
        fastq_list[name_seq[0]]=[name_seq[1], nameseq_qual[1].strip()]
    fastq_file.close()
    return fastq_list

########################################################
##Exact pattern matching, naive implementation
########################################################
from timeit import default_timer as timer

def naive(pattern, x):
    p_len = len(pattern)
    i = 0
    while i<len(x):
        if pattern == x[i:i+p_len]:
            yield i+1
        i +=1

def output_naive(fasta, fastq):
    fasta_dic = parse_fasta(fasta)
    fastq_dic = parse_fastq(fastq)
    for fq in fastq_dic:
        for fa in fasta_dic: 
            pattern = fastq_dic[fq][0]
            qual = fastq_dic[fq][1]
            x = fasta_dic[fa]
            start = timer()
            for match in naive(pattern,x):
                continue
                #print("{}\t0\t{}\t{}\t0\t{}M\t*\t0\t0\t{}\t{}".format(fq, fa,match,len(pattern),pattern,qual))
            stop = timer()
            print(pattern,len(x), stop-start)


########################################################
##Argparser
########################################################

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('fasta', help='Fasta file path for sequences')
parser.add_argument('fastq', help='Fastq file path for reads')

args = parser.parse_args()

output_naive(args.fasta,args.fastq)
