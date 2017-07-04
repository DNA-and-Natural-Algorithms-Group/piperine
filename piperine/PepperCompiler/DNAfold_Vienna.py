"""Wrapper for Vienna RNAfold."""
import os
import subprocess

import RNAfold_grammar as gram
from utils import mktemp, search_file

BREAK = "NNNNN" # Fake sequence used for strand break
def DNAfold(seq, temp, exe, par_file):
  """Runs Vienna RNAfold on sequence 'seq' at temp 'temp' (with partition function calc if 'pf' is True)."""
  infile, infilename  = mktemp(mode="w", prefix="rna_", suffix=".in")
  outfilename = infilename[:-3] + ".out"

  # Make dummy sequence which has our fake strand break (because RNAfold can't handle multistranded folding)
  dummy_seq = seq.replace("+", BREAK)

  # Write inputs
  infile.write("%s\n" % dummy_seq)
  infile.close()
  
  # Call RNAfold
  command = "%s -T %f -P %s < %s > %s" % (exe, temp, par_file, infilename, outfilename)
  subprocess.check_call(command, shell=True)
  
  # Read results
  struct, dG = gram.parseFile(outfilename)
  
  # Clean up
  os.remove(infilename)
  os.remove(outfilename)

  # Fix the structure to account for the fake strand breaking
  while BREAK in dummy_seq:
    n = dummy_seq.find(BREAK)
    dummy_seq = dummy_seq[:n] + "+" + dummy_seq[n+len(BREAK):]
    struct = struct[:n] + "+" + struct[n+len(BREAK):]
  assert seq == dummy_seq
  assert len(struct) == len(dummy_seq)
  
  return struct, dG
