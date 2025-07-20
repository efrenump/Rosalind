import numpy as np
from Bio import SeqIO

Genome = ""
Sequence = "atcaatgatcaacgtaagcttctaagcatgatcaaggtgctcacacagtttatccacaacctgagtggatgacatcaagataggtcgttgtatctccttcctctcgtactctcatgaccacggaaagatgatcaagagaggatgatttcttggccatatcgcaatgaatacttgtgacttgtgcttccaattgacatcttcagcgccatattgcgctggccaaggtgacggagcgggattacgaaagcatgatcatggctgttgttctgtttatcttgttttgactgagacttgttaggatagacggtttttcatcactgactagccaaagccttactctgcctgacatcgaccgtaaattgataatgaatttacatgcttccgcgacgatttacctcttgatcatcgatccgattgaagatcttcaattgttaattctcttgcctcgactcatagccatgatgagctcttgatcatgtttccttaaccctctattttttacggaagaatgatcaagctgctgctcttgatcatcgtttc"

for seq_record in SeqIO.parse("Vibrio cholerae DNA.fna", "fasta"):
    Genome += seq_record.seq

def PatternCount(Text, Pattern):
    Count = 0
    for i in range(0, len(Text)-len(Pattern)+1):
        if Pattern == Text[i:i+len(Pattern)]:
            Count += 1
    return Count

def FrequentWords(Text, k):
    MaxCount = 0
    FrequentPatterns = []
    for i in range(0, len(Text)-k+1):
        Pattern = Text[i:i+k]
        Count = PatternCount(Text, Pattern)
        if Count > MaxCount:
            FrequentPatterns = [Pattern]
            MaxCount = Count
        if Count == MaxCount and Pattern  not in FrequentPatterns:
            FrequentPatterns.append(Pattern)
    return FrequentPatterns

def PatternOccurrences(Genome, Pattern):
    Ocurrences = []
    for i in range(0,len(Genome)-len(Pattern)+1):
        if Genome[i:i+len(Pattern)] == Pattern:
            Ocurrences.append(i)
    return Ocurrences

def DensePatterns(Genome, L, k):
    Patterns = []
    MaxCount = 0
    for i in range(0, len(Genome)-L+1):
        FrequentPatterns = FrequentWords(Genome[i:i+L], k)
        if len(FrequentPatterns[0]) > MaxCount:
            Patterns = FrequentPatterns
            MaxCount = len(FrequentPatterns[0])
        if len(FrequentPatterns[0]) == MaxCount:
            Patterns.append(FrequentPatterns)
    return Patterns


def ReverseSequence(Sequence):
    ReverseSequence = ""
    for i in Sequence:
        if i == "a":
            ReverseSequence = "t" + ReverseSequence
        elif i == "A":
            ReverseSequence = "T" + ReverseSequence
        elif i == "t":
            ReverseSequence = "a" + ReverseSequence
        elif i == "T":
            ReverseSequence = "A" + ReverseSequence
        elif i == "c":
            ReverseSequence = "g" + ReverseSequence
        elif i == "C":
            ReverseSequence = "G" + ReverseSequence
        elif i == "g":
            ReverseSequence = "c" + ReverseSequence
        elif i == "G":
            ReverseSequence = "C" + ReverseSequence
        else:
            print(f"Sequence contains the unknown nucleotide: {i}")
    return ReverseSequence

DnaABox = FrequentWords(Sequence.upper(),9)[0]
print(DnaABox)