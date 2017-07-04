"""Wrapper for NUPACK mfe."""
import os
import subprocess

import nupack_mfe_grammar as gram
from utils import mktemp

def DNAfold(seq, temp, exe):
  """Runs NUPACK mfe on sequence 'seq' at temperature 'temp'."""
  infile, infilename  = mktemp(mode="w", prefix="mfe_", suffix=".in")
  prefix = infilename[:-3]
  outfilename = prefix+".mfe"

  seqs = seq.split("+")

  # Write inputs
  infile.write("%d\n" % len(seqs))
  for x in seqs:
    infile.write("%s\n" % x)
  for n in range(len(seqs)):
    infile.write("%d " % (n+1) )
  infile.close()
  
  # Call RNAfold
  command = "%s -T %f -material dna -multi %s" % (exe, temp, prefix)
  #print command
  subprocess.check_call(command, shell=True)
  
  # Read results
  dG, struct = gram.parseFile(outfilename)
  
  # Clean up
  os.remove(infilename)
  os.remove(outfilename)
  
  return struct, dG
