#! /usr/bin/python

desc="""Reads a set of aligned sequences in sequential FASTA format and computes several population genetics parameters"""
epilog="""Author: tgabaldon@crg.es
tgabaldon  12/08/2016
"""
#(c) tgabaldon
import string, sys, argparse, os
from collections import defaultdict
from fractions import Fraction
from math import log,sqrt


def read_fasta(inputfile): 
    #reads a fasta file and creates a dictionary
    file=open(inputfile,'r')
    lines=file.readlines()
    dictionary=defaultdict(str)
    name=''
    for line in lines:
        if line.startswith('>'):
            header = line[1:-1] # the whole header
            continue 
        dictionary[header]+=line.strip()
    return dictionary

def compareseq(seq1,seq2):
    #reads to fasta sequences and returns the number of differences and a list of  polymorphic positions
    differences=0
    polymorphic=[]
    for i in range(len(seq1)):
        if seq1[i]!=seq2[i]:
            differences=differences+1
            polymorphic.append(str(i))
    return differences,polymorphic
            
def H(n):
    """Returns an approximate value of n-th harmonic number.
       http://en.wikipedia.org/wiki/Harmonic_number
    """
    # Euler-Mascheroni constant
    gamma = 0.57721566490153286060651209008240243104215933593992
    return gamma + log(n) + 0.5/n - 1./(12*n**2) + 1./(120*n**4)      

def EH(n):
    """Returns an exact value of n-th harmonic number.
       http://en.wikipedia.org/wiki/Harmonic_number
    """
    H= lambda n: sum(Fraction(1, d) for d in xrange(1, n+1))
    return H(n)

def EHSquare(n):
    """Returns an exact value of n-th harmonic number, divided by sqaure.
       http://en.wikipedia.org/wiki/Harmonic_number
    """
    H= lambda n: sum(Fraction(1, (d*d)) for d in xrange(1, n+1))
    return H(n)

def TajimaD(Pi,S,n):
    """Returns Tajima's D
    """
    a1=float(EH(n-1))
    a2=float(EHSquare(n-1))
    n=float(n)
    b1=(n+1)/(3*(n-1))
    b2=(2*((n*n)+n+3) ) /(9*n*(n-1))
    c1=b1-(1/a1)
    c2=b2-((n+2)/(a1*n))+(a2/(a1*a1))
    e1=c1/a1
    e2=c2/((a1*a1)+a2)
    f1=Pi-(S/a1)
    f2=float(e1*S)+(e2*S*(S-1))
    #print a1,a2,b1,b2
    #print f2
    D=f1/sqrt(f2)
    return D


def main():
    parser = argparse.ArgumentParser(description=desc, epilog=epilog)

    parser.add_argument("-i", "--input", dest = "inFile", required = True,
    type = str, help = "Source file containing sequences in FASTA format")
    
    parser.add_argument("-o", "--outfile", dest = "outFile", type = str,
    default = None, help = "Output file")

    parser.add_argument("-l", "--length", dest = "seqLength", type = int,
    default = None, help = "Force sequence length. If not provided the length of the first sequence in the FASTA file will be used")

    args = parser.parse_args()

    if not os.path.isfile(args.inFile):
        sys.exit("ERROR: The input file should be defined")
    

    outfile = open(args.outFile, "w") if args.outFile else sys.stdout
    
    #reads file and creates a dictionary of the sequences
    dictionary=read_fasta(args.inFile)

    if args.seqLength:
        length = args.seqLength
    else:
        length=len(dictionary[dictionary.keys()[0]])
    # compare each pair of sequences
    global_differences=0
    global_polymorphic=[]
    compared=[]
    for header1 in dictionary:
        for header2 in dictionary:
            if header1!=header2 and header1+header2 not in compared: #this avoids comparing twice the same pair
                compared.append(header2+header1)
                seq1=dictionary[header1]
                seq2=dictionary[header2]
                differences,polymorphic=compareseq(seq1,seq2)
                #add to global counts
                global_differences=global_differences+differences
                global_polymorphic=list(set(global_polymorphic+polymorphic))
    # compute variables
    number_sequences=len(dictionary)
    mismatches=global_differences
    polymorphic_sites=len(global_polymorphic)
    comparisons=float(((number_sequences)*(number_sequences-1))/2)
    average_mismatch=global_differences/comparisons
    nucleotide_polymorphism=polymorphic_sites/float(length)
    Pi=average_mismatch/float(length)
    # uses a faster approximation if large number of sequences (set to 100)
    if number_sequences > 100:
        harmonic=float(EH(number_sequences-1))
        Theta=polymorphic_sites/harmonic
        Theta2=nucleotide_polymorphism/harmonic
    else:
        harmonic=float(H(number_sequences-1))
        Theta=polymorphic_sites/harmonic
        Theta2=nucleotide_polymorphism/harmonic


    print "Number of sequences              :", number_sequences
    print "Length                           :", length
    print "Number of segregating sites (s)  :", polymorphic_sites
    print "average number of mismatches     :", average_mismatch
    print "nucleotide polymorphism (s*)     :", nucleotide_polymorphism
    print "Nucleotide diversity (Pi)        :", Pi
    print "Theta (estimated from S)         :", Theta        
    print "Theta (estimated from S*)        :", Theta2  
    print "Theta(PI)-Theta(S)               :", average_mismatch-Theta
    print "Tajima'sD                        :", TajimaD(average_mismatch,polymorphic_sites,number_sequences)

        
if __name__ == "__main__":
    main()
