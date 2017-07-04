#!/usr/bin/env python
"""
Designs sequences randomly constrained to the forced WC-complementarity by basepairing.
Uses Joe Zadah's input and output formats for compatibility with compiler.
"""
from __future__ import division

import random
import string

from nupack_in_parser import load_design
from DNA_nupack_classes import group, complement

from ..DNAfold import DNAfold
from ..utils import error

def random_choice(group):
  """Randomly chooses an element (and gives sensical error if group has size 0)."""
  assert len(group) > 0
  return random.choice(group)

class default_list(list):
  """A list which defaults items (using a default function)."""
  def __init__(self, default_func):
    self.default = default_func
  def __getitem__(self, index):
    """Get foo[index]. If it's not in list, extend with defaults till it is."""
    if len(self) <= index:
      self.extend([self.default(a) for a in range(len(self), index+1)])
    return list.__getitem__(self, index)

class cool_set(set):
  """A set which also has some a few special attributes."""
  def __init__(self, *args):
    self.base = None
    self.const = None
    set.__init__(self, *args)

def str_replace(string, i, val):
  """Replace string[i] with val and return new string."""
  assert len(val) == 1
  return string[:i] + val + string[i+1:]

class Connect(object):
  def __init__(self):
    self.data = default_list(lambda a: default_list(lambda b: [cool_set(), cool_set([(a,b)])]))
  def add(self, struct, x, y):
    # sequence number, location on sequence, parity (answers: is not wc complement?)
    seq_x, loc_x, par_x = struct.seq_loc(x)
    seq_y, loc_y, par_y = struct.seq_loc(y)
    ### TODO: if seq_x, loc_x, par_x == seq_y, loc_y, par_y: return  # Speedup
    data_x = self.data[seq_x][loc_x]
    data_y = self.data[seq_y][loc_y]
    # Equate all involved appropriately
    data_x[par_x].update(data_y[not par_y])  
    data_x[not par_x].update(data_y[par_y])
    for i,j in data_x[True]:
      self.data[i][j][True]  = data_x[True]
      self.data[i][j][False] = data_x[False]
    for i,j in data_x[False]:
      self.data[i][j][True]  = data_x[False]
      self.data[i][j][False] = data_x[True]

def design(infilename, outfile):
  """Design sequences."""
  global d # For debugging purposes
  d = load_design(infilename)
  
  # Create complementarity matrix
  connect = Connect()
  for struct in d.structs.values():
    for x,y in struct.bonds:
      connect.add(struct, x, y)
  # Randomly color sequences
  for i, seq in enumerate(d.seqs.values()):
    for j, symb in enumerate(seq.seq):
      data = connect.data[i][j]
      if not data[True].base:
        grp = set(group[symb])
        for (x,y) in data[True]:
          grp.intersection_update(group[d.seqs.get_index(x).seq[y]])
        for (x,y) in data[False]:
          grp.intersection_update( [complement[symb] for symb in group[d.seqs.get_index(x).seq[y]]] )
        # Set constraints and choose random base
        data[True].const = tuple(grp)
        data[False].const = [complement[symb] for symb in grp]
        data[True].base = random_choice(data[True].const)
        data[False].base = complement[ data[True].base ]
      
      seq.seq = str_replace(seq.seq, j, data[True].base)

  # Test thermo and redesign to improve structure
  build_structs(d.structs)
  res = score(d.structs)
  print len(res)
  good_enough = True
  best_res = res
  while not good_enough:
    # Randomly mutate all bad structure
    for struct, index in res:
      seq_num, loc, par = struct.seq_loc(index)
      data = connect.data[seq_num][loc]
      data[True].base = random_choice(data[True].const)
      data[False].base = complement[ data[True].base ]
    # Propagate mutations
    for i, seq in enumerate(d.seqs.values()):
      for j, symb in enumerate(seq.seq):
        seq.seq = str_replace(seq.seq, j, connect.data[i][j][True].base)
    new_res = score(d.structs)
    print len(new_res), len(res), len(best_res)
    res = new_res
    best_res = min(best_res, res)
    #raw_input()
    
  output_sequences(d, connect, outfile)

def build_structs(structs):
  for struct in structs.values():
    breaks = [i for i, symb in enumerate(struct.struct) if symb == "+"]
    lengths = [y - x - 1 for x, y in zip([-1] + breaks, breaks + [len(struct.struct)])] # = diffs of breaks
    #print breaks, lengths
    lengths = iter(lengths)
    # Build strands list
    strand = []
    struct.strands = [strand]
    l = lengths.next()
    #print struct.seqs
    for seq in struct.seqs:
      if l == 0:
        l = lengths.next()
        strand = []
        struct.strands.append(strand)
      if seq.length <= l:
        strand.append(seq)
        #print struct.name, seq.name, l, type(l)
        l -= seq.length
      else:
        #print struct.name, seq.name, l, seq.length
        raise ValueError, "Sequences do not match structure"
    assert l == 0 and end_iter(lengths)

def end_iter(foo):
  """Iterate to the end of foo."""
  try:
    foo.next()
    return False
  except StopIteration:
    return True

def build_seq(struct):
  """Produce flat sequence for a structure by expanding out the strands and subsequences/subdomains."""
  return string.join([ string.join([seq.seq for seq in strand], "") for strand in struct.strands], "+")

def score(structs):
  """Find unwanted base-pairing"""
  diffs = []
  for struct in structs.values():
    struct.seq = build_seq(struct)
    #print struct.name, struct.seq
    struct.mfe_struct, dG = DNAfold(struct.seq)
    diffs += diff(struct, struct.struct, struct.mfe_struct)
  return diffs

def diff(lable, str_a, str_b):
  """Return bp differences between 2 secondary structures."""
  ### TODO: make correct by having:
  ###       ..(((...+...)))..  to  ..(((..)+(..)))..  be 4 errors not just 2
  #print lable.name, str_a, str_b
  assert len(str_a) == len(str_b)
  str_a = str_a.replace("+", "")
  str_b = str_b.replace("+", "")
  return [(lable, i) for i, (a, b) in enumerate(zip(str_a, str_b)) if a != b]

def output_sequences(d, connect, fn):
  """Output sequences in Joe's format."""
  f = open(fn, "w")
  for struct in d.structs.values():
    # Write structure (with dummy content)
    f.write("%d:%s\n" % (0, struct.name))
    gc_content = (struct.seq.count("C") + struct.seq.count("G")) / (len(struct.seq) - len(struct.strands)) # - len(struct.strands) because of the +'s in the struct.seq
    f.write("%s %f %f %d\n" % (struct.seq, 0, gc_content, 0))
    f.write("%s\n" % struct.struct)   # Target structure
    f.write("%s\n" % struct.mfe_struct)  # MFE structure of chosen sequences
  for seq in d.seqs.values():
    # Write sequence (with dummy content)
    f.write("%d:%s\n" % (0, seq.name))
    gc_content = (seq.seq.count("C") + seq.seq.count("G")) / seq.length
    f.write("%s %f %f %d\n" % (seq.seq, 0, gc_content, 0))
    f.write(("."*seq.length+"\n")*2) # Write out dummy structures.
    # Write wc of sequence (with dummy content)
    seq = seq.wc
    f.write("%d:%s\n" % (0, seq.name))
    #wc_seq = string.join([complement[symb] for symb in seq.seq[::-1]], "")
    f.write("%s %f %f %d\n" % (seq.seq, 0, gc_content, 0))
    f.write(("."*seq.length+"\n")*2) # Write out dummy structures.
  f.write("Total n(s*) = %f" % 0)
  f.close()

if __name__ == "__main__":
  import sys
  import re
  
  from find_file import find_file
  
  if sys.version_info < (2, 5):
    error("Must use python 2.5 or greater.")
    
  infilename = find_file(sys.argv[1])
  
  # Infer the basename if a full filename is given
  basename = infilename
  p = re.match(r"(.*)\.des\Z", basename)
  if p:
    basename = p.group(1)
  # Set output file name
  outfilename = basename + ".mfe"
  
  design(infilename, outfilename)
