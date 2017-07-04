#!/usr/bin/env python

# Reports the longest exact matches in a sequence.
# Intended for use to examine what are the worst offenders that spuriousSSM is reporting.
# Typical use: Copy designed sequence to a new file.  Invoke findspurious with that filename.
# Example: ./findspurious examples/TAE.sp

# Note that, unlike spuriousSSM, the "intra-complex" category here is not necessarily "inter-strand",
# and that the most general category here is not necessarily either "inter-strand" nor "inter-complex".
# I.e. here each category considers a superset of the more restricted category, rather than being mutually exclusive as in spuriousSSM.

# Erik Winfree, Feb 2015
 
import re
import sys

if len(sys.argv) != 2:
    sys.exit("Usage:  findspurious <sequencefilename> \n")

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

def WC(s):
    return s[::-1].replace('A','t').replace('T','a').replace('C','g').replace('G','c').replace('a','A').replace('t','T').replace('c','C').replace('g','G')

k=0
found=1
while found:
    k=k+1
    found=0
    for i in range(1,len(lines)-k):
        word = lines[i:i+k]
        if dna.match(word).group()==word:
            if lines[i+1:].find(word) != -1:
                matchword=word
                found=1
k=k-1
print "Longest exact match: " + matchword + " of length " + str(len(matchword))

# this detects self-complementary words
k=0
found=1
while found:
    k=k+1
    found=0
    for i in range(1,len(lines)-k):
        word = WC(lines[i:i+k])
        if dna.match(word).group()==word:
            if lines[i:].find(word) != -1:
                matchword=word
                found=1
k=k-1
print "Longest complementary match: " + matchword + " and " + WC(matchword)  + " of length " + str(len(matchword))

k=0
found=1
while found:
    k=k+1
    found=0
    for i in range(1,len(lines)-k):
        word = lines[i:i+k]
        if dna.match(word).group()==word:
            strandend1=lines[i+1:].find(' ')
            strandend2=lines[i+1:].find('\n')
            if strandend1==-1 and strandend2==-1:
                rest=lines[i+1:]
            else:
                rest=lines[i+1:max(strandend1,strandend2)]
            if rest.find(word) != -1:
                matchword=word
                found=1
k=k-1
print "Longest intra-strand exact match: " + matchword  + " of length " + str(len(matchword))

# this detects self-complementary words
k=0
found=1
while found:
    k=k+1
    found=0
    for i in range(1,len(lines)-k):
        word = WC(lines[i:i+k])
        if dna.match(word).group()==word:
            strandend1=lines[i+1:].find(' ')
            strandend2=lines[i+1:].find('\n')
            if strandend1==-1 and strandend2==-1:
                rest=lines[i:]
            else:
                rest=lines[i:max(strandend1,strandend2)]
            if rest.find(word) != -1:
                matchword=word
                found=1
k=k-1
print "Longest intra-strand complementary match: " + matchword + " and " + WC(matchword)  + " of length " + str(len(matchword))



k=0
found=1
while found:
    k=k+1
    found=0
    for i in range(1,len(lines)-k):
        word = lines[i:i+k]
        if dna.match(word).group()==word:
            strandend1=lines[i+1:].find('  ')
            strandend2=lines[i+1:].find('\n\n')
            if strandend1==-1 and strandend2==-1:
                rest=lines[i+1:]
            else:
                rest=lines[i+1:max(strandend1,strandend2)]
            if rest.find(word) != -1:
                matchword=word
                found=1
k=k-1
print "Longest intra-complex exact match: " + matchword  + " of length " + str(len(matchword))

# this detects self-complementary words
k=0
found=1
while found:
    k=k+1
    found=0
    for i in range(1,len(lines)-k):
        word = WC(lines[i:i+k])
        if dna.match(word).group()==word:
            strandend1=lines[i+1:].find('  ')
            strandend2=lines[i+1:].find('\n\n')
            if strandend1==-1 and strandend2==-1:
                rest=lines[i:]
            else:
                rest=lines[i:max(strandend1,strandend2)]
            if rest.find(word) != -1:
                matchword=word
                found=1
k=k-1
print "Longest intra-complex complementary match: " + matchword + " and " + WC(matchword)  + " of length " + str(len(matchword))

print "\n"
