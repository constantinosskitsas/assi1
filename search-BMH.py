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
##Exact pattern matching, BMH
########################################################

def BMH(x,p):
    m = len(p)
    #build jump table with unknown alphabet
    jump = {}
    for idx in range(m-1):
        jump[p[idx]] = idx
    for i in x:
        if i not in jump:
            jump[i] = -1
    #find and report matches
    j=m-1
    i=m-1
    while j<len(x):
        if x[j-(m-1)+i] == p[i]:
            if i == 0:
                yield (j-(m-1)+1)
                j = j+(m-1)-jump[x[j]]
                i=m-1
            else:
                i-=1
        else:
            j = j+(m-1)-jump[x[j]]
            i = m-1

def output_BMH(fasta,fastq):
    fasta_dic = parse_fasta(fasta)
    fastq_dic = parse_fastq(fastq)
    for fq in fastq_dic:
        for fa in fasta_dic:
            pattern = fastq_dic[fq][0]
            qual = fastq_dic[fq][1]
            x = fasta_dic[fa]
            for match in BMH(x,pattern):
                print("{}\t0\t{}\t{}\t0\t{}M\t*\t0\t0\t{}\t{}".format(fq, fa,match,len(pattern),pattern,qual))


########################################################
##Argparser
########################################################

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('fasta', help='Fasta file path for sequences')
parser.add_argument('fastq', help='Fastq file path for reads')

args = parser.parse_args()

output_BMH(args.fasta,args.fastq)