#!/usr/bin/python3 -u
import math
import sys
from os.path import basename
from statistics import variance

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import GC, gc_fraction

arguments = sys.argv[1:]

print("File\tA percent\tT percent\tC percent\tG percent\tGC percent\tAT/GC ratio\tNucleScore\tATG\tTGA\tTAG\tTAA\tGenome size")

for argv in arguments:
    #print("File: "+argv+"\n\n")
    file = argv
    globalSeq = ''
    fasta_sequences = SeqIO.parse(open(file), 'fasta')

    for seq_record in fasta_sequences:
        globalSeq += str(seq_record.seq)

    gcpercent = gc_fraction(globalSeq) * 100

    my_dna = Seq(globalSeq)

    ade = my_dna.count("A")
    thy = my_dna.count("T")
    gua = my_dna.count("G")
    cyt = my_dna.count("C")
    n = my_dna.count("N")
    length = len(globalSeq)

    aPercent = (ade/length)*100
    tPercent = (thy/length)*100
    gPercent = (gua/length)*100
    cPercent = (cyt/length)*100
    nPercent = (n/length) * 100

    atgcRatio = (ade + thy) / (gua + cyt)
    percentList = (aPercent, tPercent, gPercent, cPercent, nPercent)
    variance_value = variance(percentList)
    nucleScore = math.log2((variance_value * gcpercent * atgcRatio ** 3) / math.sqrt(length))

    # act = my_dna.find("ACT")

    atg = my_dna.count('ATG')
    tga = my_dna.count('TGA')
    tag = my_dna.count('TAG')
    taa = my_dna.count('TAA')

    label = basename(file)

    # Summary file
    print(label + "\t" + str(
            aPercent) + "\t" + str(tPercent) + "\t" + str(cPercent) + "\t" + str(gPercent) + "\t" + str(
            gcpercent) + "\t" + str(atgcRatio) + "\t" + str(nucleScore) + "\t" + str(atg) + "\t" + str(tga) + "\t" + str(tag) + "\t" + str(taa) + "\t" + str(length))
