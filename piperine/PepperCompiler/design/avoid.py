"""
The k-sequence avoidance algorithm attempts to construct sequences which avoid 
having any subsequences of length k which are complementary unless explicitly
forced to be so.

Considered extension:
  Allow complementarity on overlapping regions.
  Thus k=3 would allow ACGT even though ACG ~ CGT

Uses the input format from Winfree's SpuriousC algorithm.
"""

import sys
import random

import DNA_nupack_classes as DNA_classes

def last(n, foo):
  """Gets the last n items in foo."""
  if n == 0:
    return foo[len(foo):]
  else:
    return foo[-n:]  # Note that foo[-0:] would fail.

def get_group(s):
  """Convert a letter into a group. i.e. "S" -> "CG", "N" -> "ATCG", etc. """
  if s == " ":
    return " "
  else:
    return DNA_classes.group[s]

class Design(object):
  def __init__(self, st, eq, wc):
    """Load in the parameters for the system."""
    assert len(st) == len(eq) == len(wc)
    self.n = len(st)
    self.st = map(get_group, st)  # Get the groups from the letter.
    self.eq = eq
    self.wc = wc
  
  def avoid(self, k):
    """
    Avoid k-subsequence repeats or complementarity 
    with constraints that nt be in st, and 
    eq and wc specify representatives for equality and complementarity.
    """
    # We will recurse at least self.n times.
    # This might fail on to.dna for >25,000 recursion with Seg Fault!
    sys.setrecursionlimit(max(10*self.n, 1000))
    return self.avoid_rec(k, 0, "", {})
    
  def avoid_rec(self, k, i, part_seq, bad_seqs):
    """Recursively steps through the algorithm at position i avoiding bad_sequences."""
    # If we've assigned all bases, we're done.
    if i >= self.n:
      return part_seq
    
    # If it's a strand break, pop it on and keep going.
    if self.st[i] == " ":
      return self.avoid_rec(k, i+1, part_seq + " ", bad_seqs)
    
    # If we've already assigned an equal or wc constraint we must respect it
    if self.eq[i] < i:
      order = [part_seq[self.eq[i]]] # There is only one choice
    elif self.wc[i] != None and self.wc[i] < i:
      order = [DNA_classes.complement[part_seq[self.wc[i]]]] # There is only one choice
    
    # Otherwise, randomly order the allowed nts
    else:
      order = list(self.st[i])
      random.shuffle(order)

    # Try each nucleotide in group until one works
    for nt in order:
      new_end = last((k-1), part_seq) + nt  # End of new seq with nucleotide added
      # If the lask k nts are in a single strand
      if " " not in new_end:
        # If the last 'k' nts are eq to a past seq ...
        if new_end in bad_seqs:
          j = bad_seqs[new_end]
          # ... and they aren't suposed to be, fail
          if self.eq[i+1-k:i+1] != self.eq[j+1-k:j+1]:
            #print "eq", j, i, part_seq + nt
            continue
        
        comp_end = DNA_classes.seq_comp(new_end)  # Complement
        # If the last 'k' nts are wc to a past seq ...
        if comp_end in bad_seqs:
          j = bad_seqs[comp_end]
          # ... and they aren't suposed to be, fail.
          # Note: we have reversed the indexing to wc.
          # TODO-test: ignore overlapping regions (j > i-k) because they will not bond.
          #if j <= i-k and self.wc[i:i-k:-1] != self.eq[j+1-k:j+1]:
          if self.wc[i:i-k:-1] != self.eq[j+1-k:j+1]:
            #print "wc", j, i, part_seq + nt
            continue
      
      # Otherwise, use nt and step deeper.
      #new_bad = bad_seqs.copy() # Don't mutate bad_seqs
      
      # We mutate bad_seqs and then correct it.
      modify = new_end not in bad_seqs
      if modify:
        bad_seqs[new_end] = i
      
      # Try using this nt
      res = self.avoid_rec(k, i+1, part_seq + nt, bad_seqs)
      
      # Correct bad_seqs
      if modify:
        del bad_seqs[new_end]
      
      if res:
        return res
      # Otherwise continue trying
    
    return None # All nucleotides fail
  

def testU(k, n):
  """Test for a length n unpaired single strand."""
  d = Design("N"*n, range(n), [-1]*n)
  return d.avoid(k)

def testH(k, n):
  """Test for a length n helix."""
  d = Design("N"*n + " " + "N"*n, range(2*n+1), range(2*n, -1, -1))
  return d.avoid(k)
