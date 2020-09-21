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
##Exact pattern matching, KMP
########################################################

def border_array_lin(x):
    B = [0]*len(x)
    bax = [0]*len(x)
    #make a border array
    for i in range(1,len(x)):
        b = B[i-1]
        while b>0 and x[i] != x[b]:
            b = B[b-1]
        if x[i] == x[b]:
            B[i] = b+1
            b+=1
        else:
            B[i] = 0
            b=0
        #make bax - longest border where the following characters differ
        if i == len(x)-1:
            bax[i] = b
        elif x[i+1] != x[b]:
            bax[i] = b
        else:
            bax[i] = bax[B[b-1]]

    return bax

def KMP(x, p):
    bax = border_array_lin(p)
    m = len(p)-1
    i=0
    j=0
    while j<len(x)-m+i:
        if x[j] == p[i]:
            if i == m:
                yield (j-m+1)
                j+=1
                i = bax[i]
            else:
                j+=1
                i+=1
        elif i == 0:
            j+=1
        else:
            i = bax[i-1]

def output_KMP(fasta,fastq):
    fasta_dic = parse_fasta(fasta)
    fastq_dic = parse_fastq(fastq)
    for fq in fastq_dic:
        for fa in fasta_dic:
            pattern = fastq_dic[fq][0]
            qual = fastq_dic[fq][1]
            x = fasta_dic[fa]
            for match in KMP(x,pattern):
                print("{}\t0\t{}\t{}\t0\t{}M\t*\t0\t0\t{}\t{}".format(fq, fa,match,len(pattern),pattern,qual))



########################################################
##Argparser
########################################################

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('fasta', help='Fasta file path for sequences')
parser.add_argument('fastq', help='Fastq file path for reads')

args = parser.parse_args()

output_KMP(args.fasta,args.fastq)