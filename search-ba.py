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
##Exact pattern matching, border array
########################################################

def border_array(pattern,x):
    pattern1 = pattern +"$"
    #border array for pattern 
    BP = [0]*len(pattern1)
    for i in range(1,len(pattern1)):
        b = BP[i-1]
        while b>0 and pattern1[i] != pattern1[b]:
            b = BP[b-1]
        if pattern1[i] == pattern1[b]:
            BP[i] = b+1
        else:
            BP[i] = 0
    #border array for pattern$x
    B = [0]*len(x)
    for i in range(0,len(x)):
        if i == 0:
            b=0
        else:
            b = B[i-1]
        while b>0 and x[i] != pattern1[b]:
            b = BP[b-1]
        if x[i] == pattern1[b]:
            B[i] = b+1
        else:
            B[i] = 0
    for i in range(len(B)):
        if B[i] == len(pattern):
            yield (i-(len(pattern)-1)+1) 

def output_borderarray(fasta, fastq):
    fasta_dic = parse_fasta(fasta)
    fastq_dic = parse_fastq(fastq)
    for fq in fastq_dic:
        for fa in fasta_dic:
            pattern = fastq_dic[fq][0]
            qual = fastq_dic[fq][1]
            x = fasta_dic[fa]
            for match in border_array(pattern,x):
                print("{}\t0\t{}\t{}\t0\t{}M\t*\t0\t0\t{}\t{}".format(fq, fa,match,len(pattern),pattern,qual))


########################################################
##Argparser
########################################################

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('fasta', help='Fasta file path for sequences')
parser.add_argument('fastq', help='Fastq file path for reads')

args = parser.parse_args()

output_borderarray(args.fasta,args.fastq)