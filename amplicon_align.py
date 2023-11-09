#!/usr/bin/env python3
import magnumopus
import argparse

parser = argparse.ArgumentParser(description = "Perform in-silico PCR on two assemblies and align the amplicons")
parser.add_argument('-1', metavar='ASSEMBLY1', type=str, help='Path to first assembly file')
parser.add_argument('-2', metavar='ASSEMBLY2', type=str, help='Path to second assembly file')
parser.add_argument('-p', metavar='PRIMERS', type=str, help='Path to the primer file')
parser.add_argument('-m', metavar='MAX_AMPLICON_SIZE', type=str, help='maximum amplicon size for isPCR')
parser.add_argument('--match', metavar='MATCH', type=str, help='match score to use in alignment')
parser.add_argument('--mismatch', metavar='MISMATCH', type=str, help='mismatch penalty to use in alignment')
parser.add_argument('--gap', metavar='GAP', type=str, help='gap penalty to use in alignment')
args = parser.parse_args()

#Get the inputs as variables
Assembly1 = args.__dict__['1']
Assembly2 = args.__dict__['2']
Primers = args.__dict__['p']
Max_Amplicon_size = args.__dict__['m']
Match = args.__dict__['match']
Mismatch = args.__dict__['mismatch']
Gap = args.__dict__['gap']

amplicons1 = magnumopus.ispcr(Primers, Assembly1, int(Max_Amplicon_size))
amplicons2 = magnumopus.ispcr(Primers, Assembly2, int(Max_Amplicon_size))

#Getting our amplicons
amplicons1 = (amplicons1.split("\n"))[1]
amplicons2 = (amplicons2.split("\n"))[1]

#Checking reverse direction
complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
reverse_complement_1 = "".join(complement.get(base, base) for base in reversed(amplicons2))

#Do nw algorithm
aln1, score1 = magnumopus.needleman_wunsch(amplicons1, amplicons2, int(Match), int(Mismatch), int(Gap))
aln2, score2 = magnumopus.needleman_wunsch(amplicons1, reverse_complement_1, int(Match), int(Mismatch), int(Gap))

#Print best match
if score1 > score2:
    print("\n".join(aln1))
    print(score1)
else:
    print("\n".join(aln2))
    print(score2)   



