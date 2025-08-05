import numpy as np
from Bio import SeqIO
import matplotlib.pyplot as plt
import time

Genome = ""

Nucleotides = "ACGT"
NucleotidesInt = "0123"
NucleotideTable = str.maketrans(Nucleotides, NucleotidesInt)

for seq_record in SeqIO.parse("Escherichia coli DNA.fasta", "fasta"):
    Genome += str(seq_record.seq)

def PatternCount(Text, Pattern):
    Count = Start = 0
    while True:
        Start = Text.find(Pattern, Start) + 1
        if Start > 0:
            Count += 1
        else:
            return Count

def PatternOccurrences(Text, Pattern):
    Ocurrences = []
    Start = 0

    while True:
        Start = Text.find(Pattern, Start) + 1
        if Start > 0:
            Ocurrences.append(Start-1)
        else:
            return Ocurrences

def PatternToNumber(Pattern):
    return int(Pattern.translate(NucleotideTable), 4)

def NumberToPattern(Number, k=-1):
    if k==-1:
        k=int(np.floor(np.log10(Number)/np.log10(4))+1)

    if Number < 0 or Number >= 4**k:
        return 0
    else:
        out = ["A"] * k
        i = k-1

        while i>=0: #Euclid's algorithm but with bitwise operators
            out[i] = Nucleotides[Number & 3] #Same as taking mod 4.
            Number >>= 2 #Same as removing the last character.
            i -= 1

        return "".join(out)

def ComputeFrequencies(Text, k): #O(n)
    Frequencies = [0] * 4**k
    Mask = len(Frequencies) - 1

    Pattern = PatternToNumber(Text[0:k])
    Frequencies[Pattern] += 1

    for i in range(k,len(Text)):
        Pattern = ((Pattern << 2) | PatternToNumber(Text[i])) & Mask  
        Frequencies[Pattern] += 1  

    return Frequencies

def FrequentWords(Text, k): #O(n + 4^k)
    Frequencies = ComputeFrequencies(Text,k)
    MaxCount = max(Frequencies)

    return [NumberToPattern(i,k) for i, v in enumerate(Frequencies) if v == MaxCount]

def EncodePatterns(Text, k):
    Pattern = PatternToNumber(Text[0:k])
    EncodedPatterns = [Pattern]
    Mask = 4**k - 1

    for i in range(k, len(Text)):
        Pattern = ((Pattern << 2) | PatternToNumber(Text[i])) & Mask
        EncodedPatterns.append(Pattern)

    return EncodedPatterns

def FrequentWords2(Text, k): #O(nlog(n)), better for k>=13
    Index = EncodePatterns(Text, k)
    Count = [1]*(len(Text)-k+1)

    Index.sort()
    for i in range(1,len(Text)-k+1):
        if Index[i] == Index[i-1]:
            Count[i] = Count[i-1] + 1

    MaxCount = max(Count)
    return [NumberToPattern(Index[i], k) for i in range(0,len(Text)-k+1) if Count[i] == MaxCount]

def DenseWords(Text, k, t, L):
    Frequencies = ComputeFrequencies(Text[0:L], k)
    Clump = [(Frequencies[i] >= t) for i in range(0, 4**k)]

    for i in range(0, len(Text)-L+1):
        Frequencies[PatternToNumber(Text[i:i+k])] -= 1
        LastPattern = PatternToNumber(Text[i+L-k+1:i+L+1])
        Frequencies[LastPattern] += 1
        if Frequencies[LastPattern] >= t:
            Clump[LastPattern] = 1

    return [NumberToPattern(i,k) for i in range(0, 4**k) if Clump[i] == 1]

CompNucleotideTable = str.maketrans(Nucleotides, Nucleotides[::-1])

def ReverseSequenceStr(Sequence):
    ReverseSequence = Sequence.translate(CompNucleotideTable)
    
    return ReverseSequence[::-1]

CompNucleotideTableInt = str.maketrans(NucleotidesInt, NucleotidesInt[::-1])

def ReverseSequenceInt(Sequence, k):
    Result = 0
    for i in range(0, k):
        Result = (Result << 2) | ((Sequence & 3) ^ 3)
        Sequence >>= 2
    return Result

def SkewGC(Genome, i): #Increases in the forward (lagging) half-strand and decreases in the reverse (leading) half-strand
    return PatternCount(Genome[0:i].upper(), "G") - PatternCount(Genome[0:i].upper(), "C")

def SkewGCArray(Genome):
    SkewGC = 0
    SkewValues = []
    for i in range(0, len(Genome)):
        if Genome[i] == "G": 
            SkewGC += 1
        elif Genome[i] == "C":
            SkewGC -= 1
        SkewValues.append(SkewGC)
    return SkewValues

def HammingDistanceStr(Word1, Word2):
    if len(Word1) == len(Word2):
        return sum(c1 != c2 for c1, c2 in zip(Word1, Word2))
    else:
        print("Both sequences don't have the same length")
        return -1

def HammingDistanceInt(Word1, Word2, k):
    Result = Word1 ^ Word2
    Result |= Result >> 1
    Result &= (2 ** (2*k)) // 3

    return Result.bit_count()

def NeighborsStr(Pattern: str, d):
    Neighborhood = []

    if d == 0:
        return [Pattern]
    elif len(Pattern) == 1:
            return ["A", "C", "G", "T"]
    else:
        Suffix = Pattern[1:len(Pattern)]
        SuffixNeighborhood = NeighborsStr(Suffix, d)

        for i in SuffixNeighborhood:
            if HammingDistanceStr(Suffix, i) < d:
                for j in Nucleotides:
                    Neighborhood.append(j + i)
            else:
                Neighborhood.append(Pattern[0] + i)

    return Neighborhood

def NeighborsInt(Pattern: int, k, d):
    Neighborhood = []

    if d == 0:
        return [Pattern]
    elif k == 1:
        return [0, 1, 2, 3]
    else:
        Mask = 4**(k-1)-1
        Suffix = Pattern & Mask
        SuffixNeighborhood = NeighborsInt(Suffix, k-1, d)

        for i in SuffixNeighborhood:
            if HammingDistanceInt(Suffix, i, k-1) < d:
                for j in [0, 1, 2, 3]:
                    Neighborhood.append((j << (2*(k-1))) | i)
            else:
                Neighborhood.append((Pattern & ~Mask) | i)

    return Neighborhood

def ApproximatePatternCount(Text, Pattern, d): ##Work with numbers instead of strings
    Count = 0
    for i in range(0, len(Text)-len(Pattern)+1):
        Window =  Text[i:i+len(Pattern)]
        HammingDistance = 0
        for c1, c2 in zip(Pattern, Window):
            if c1 != c2:
                HammingDistance += 1
            if HammingDistance > d:
                break
        if HammingDistance <= d:
            Count += 1
                
    return Count

def ComputeApproximateFrequencies(Text, k, d):
    ApproximateFrequencies = [0] * 4**k

    Frequencies = [0] * 4**k

    Mask = len(Frequencies) - 1

    Pattern = PatternToNumber(Text[0:k])
    Frequencies[Pattern] += 1
    Seen = []
    SeenFlag = [0] * 4**k

    SeenFlag[Pattern] = 1
    Seen.append(Pattern)

    for i in range(k,len(Text)):
        Pattern = ((Pattern << 2) | PatternToNumber(Text[i])) & Mask  
        Frequencies[Pattern] += 1
        if SeenFlag[Pattern] == 0:
            SeenFlag[Pattern] = 1
            Seen.append(Pattern)

    for i in Seen:
        for j in NeighborsInt(i, k, d):
            ApproximateFrequencies[j] += Frequencies[i]

    return ApproximateFrequencies

def FrequentApproximateWords(Text, k, d):
    ApproximateFrequencies = [0] * 4**k

    Frequencies = [0] * 4**k

    Mask = len(Frequencies) - 1

    Pattern = PatternToNumber(Text[0:k])
    Frequencies[Pattern] += 1
    Seen = []
    SeenFlag = [0] * 4**k

    SeenFlag[Pattern] = 1
    Seen.append(Pattern)

    for i in range(k,len(Text)):
        Pattern = ((Pattern << 2) | PatternToNumber(Text[i])) & Mask  
        Frequencies[Pattern] += 1
        if SeenFlag[Pattern] == 0:
            SeenFlag[Pattern] = 1
            Seen.append(Pattern)
     
    Max = 0
    MostFrequent = []

    for i in Seen:
        for j in NeighborsInt(i, k, d):
            ApproximateFrequencies[j] += Frequencies[i]

            NeighborCount = ApproximateFrequencies[j]
            if NeighborCount > Max:
                MostFrequent = [j]
                Max = NeighborCount
            elif NeighborCount == Max:
                if j not in MostFrequent:
                    MostFrequent.append(j)

    return [NumberToPattern(i,k) for i in MostFrequent]

def FrequentApproximateComplementaryWords(Text, k, d):
    Frequencies = ComputeApproximateFrequencies(Text, k, d)
    ComplementaryFrequencies = [0] * 4**k

    Max = 0
    MostFrequent = []

    for i in range(0, 4**k):
        Reverse = ReverseSequenceInt(i, k)
        ComplementaryFrequencies[i] = Frequencies[i]
        if Reverse != i:
            ComplementaryFrequencies[i] += Frequencies[Reverse]

        Count = ComplementaryFrequencies[i]
        if Count > Max:
            Max = Count
            MostFrequent = [i]
        elif Count == Max:
            MostFrequent.append(i)

    return [NumberToPattern(i, k) for i in MostFrequent]

SkewArray = SkewGCArray(Genome)
MinSkewY = min(SkewArray)
MinSkewX = [i for i, v in enumerate(SkewArray) if v == MinSkewY]
print(MinSkewX)

fig, ax = plt.subplots()
ax.plot(SkewArray)
ax.vlines(MinSkewX[0], 0, 1, transform=ax.get_xaxis_transform(), colors="red")

start = time.time()
print(FrequentApproximateComplementaryWords(Genome[MinSkewX[0]-100:MinSkewX[0]+400], 9, 1))
end = time.time()
print(end - start)

plt.show()

