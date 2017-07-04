"""Multistrand wrapper"""
from __future__ import division

import os
import subprocess
import re

OUTFILE = "Logfile"

def DNAkinfold(base_infile, params):
  """
  Expand (or replace) parameters in base_infile with those supplied from 
  the command line.
  """
  in_name  = base_infile + ".in"
  
  fbase = open(base_infile, "r")
  fin = open(in_name, "w")
  
  # Copy all the file content that's not single line parameters
  #   and save all single line parameters to the parameter list (or ignore).
  for line in fbase:
    p = re.match(r"#([^=]*)=(.*)\n", line)
    if p: # If it's a single line parameter
      name, value = p.group(1, 2)
      if name not in params:  # We ignore overridden parameters.
        params[name] = value
    else:
      fin.write(line)
  
  fbase.close()
  
  if OUTFILE not in params:
    params[OUTFILE] = base_infile + ".out"
  out_name = params[OUTFILE]
  
  # Write out all parameters.
  for name, value in params.items():
    fin.write("#%s=%s\n" % (name, value))
  
  fin.close()
  
  if os.path.isfile(out_name):
    os.remove(out_name)
  # Run Multistrand!
  command = "Multistrand < %s" % in_name
  if params["OutputInterval"] == -1: # If we say quiet, we mean it.
    command += " > /dev/null"
    
  print command
  subprocess.check_call(command, shell=True)

if __name__ == "__main__":
  import sys
  
  try:
    infile = sys.argv[1]
    params = dict([arg.split("=", 1) for arg in sys.argv[2:]])
  except:
    print "usage: python multihelp.py infile [ARG=VAL ...]"
    sys.exit(1)
  DNAkinfold(infile, params)
