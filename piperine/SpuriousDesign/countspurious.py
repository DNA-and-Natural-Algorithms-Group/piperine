#!/usr/bin/env python

import re
import sys

if len(sys.argv) != 2:
    sys.exit("Usage:  countspurious.py <sequencefilename> \n")

print "\n"

filename = sys.argv[1]
text_file = open(filename, "r")
# lines = text_file.readlines()
lines = text_file.read()
print "Scanning sequences:\n"+lines
# print len(lines)
print "\n"
text_file.close()

dna=re.compile('[ACTG]*')

# returns the WC string
def WC(s):
    return s[::-1].replace('A','a').replace('T','t').replace('C','c').replace('G','g').replace('a','T').replace('t','A').replace('c','G').replace('g','C')

kmax=20
spurious=[0]*kmax
spuriousWC=[0]*kmax

for k in range(1,kmax+1):
   subseq=set()
   count=0
   countWC=0
   for i in range(1,len(lines)):
       word = lines[i:i+k]
       if dna.match(word).group()==word:
           if word in subseq:
               count+=1
           subseq.add(word)
           if WC(word) in subseq:
               countWC+=1
   spurious[k-1] = count
   spuriousWC[k-1] = countWC

print "Counts of exact matches of length 1 through " + str(kmax) + " : " + str(spurious)
print "Counts of complementary matches of length 1 through " + str(kmax) + " : " + str(spuriousWC)
