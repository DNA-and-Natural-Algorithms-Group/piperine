#!/usr/bin/env python
"""
Designs sequences using a k-sequence avoiding algorithm.
Uses PIL input and Zadeh's .des output formats for compatibility with compiler.
"""

import avoid
from constraint_load import Convert

def design(in_name, k):
  convert = Convert(in_name, struct_orient=False)
  eq, wc, st = convert.get_constraints()
  
  def st_map(x):
    return x  if x != None  else " "
  st = map(st_map, st)
  
  d = avoid.Design(st, eq, wc)
  print d.avoid(k)
  # TODO: something with these sequences.

if __name__ == "__main__":
  import sys
  
  try:
    in_name = sys.argv[1]
    k = int(sys.argv[2])
  except:
    print "Usage: python avoid_design.py infilename k"
    sys.exit(1)
  
  design(in_name, k)
