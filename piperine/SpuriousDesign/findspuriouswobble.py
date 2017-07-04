#!/usr/bin/env python

import re
import sys

if len(sys.argv) != 2:
    sys.exit("Usage:  findspuriouswobble <sequencefilename> \n")

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

# returns the WC string, no wobble.
def WC(s):
    return s[::-1].replace('A','a').replace('T','t').replace('C','c').replace('G','g').replace('a','T').replace('t','A').replace('c','G').replace('g','C')


# creates a regular expression that matches WC-with-wobbles strings
def WCwobble(s):
    return s[::-1].replace('A','a').replace('T','t').replace('C','c').replace('G','g').replace('a','T').replace('t','[AG]').replace('c','G').replace('g','[CU]')


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


# this detects self-complementary subsequences
k=0
match=True
while match:
    k=k+1
    pat=WC(lines[1:1+k])
    for i in range(2,len(lines)-k):
        word = WC(lines[i:i+k])
        pat = pat+'|'+word
    p=re.compile(pat)
    match=p.search(lines)
    if match:
        matchword=match.group()
k=k-1
print "Longest complementary match: " + matchword + " of length " + str(len(matchword))

# this detects self-complementary subsequences
k=0
match=True
while match:
    k=k+1
    pat=WCwobble(lines[1:1+k])
    for i in range(2,len(lines)-k):
        word = WCwobble(lines[i:i+k])
        pat = pat+'|'+word
    p=re.compile(pat)
    match=p.search(lines)
    if match:
        matchword=match.group()
k=k-1
print "Longest wobble-complementary match: " + matchword + " of length " + str(len(matchword))
