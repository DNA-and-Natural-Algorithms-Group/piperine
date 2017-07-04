#!/usr/bin/env python2.6

import sys
from utils import error
if sys.version_info < (2, 5):
  error("Must use python 2.5 or greater.")

import re
import os
import pickle
import time

from system_class import load_file
from utils import match, warning, error

def parse_fixed(line):
  """Parse a line in the fixed file."""
  m = match(r"(\w+) ([\w_-]+)[ \t]*=[ \t]*([ATCGNS]+)(?: #.*)?", line)
  try:
      type_, name, seq = m.groups()
  except AttributeError:
      raise ValueError(line)
  return type_, name, seq

def load_fixed(filename):
  """Load a file of sequences to fix."""
  if not os.path.isfile(filename):
    error("Cannot fix sequences. No such file '%s'." % filename)
  
  f = open(filename, "r")
  return [parse_fixed(line) for line in f if not re.match(r"\s*(#.*)?\s*\Z", line)]
  
def compiler(basename, args, outputname, savename, fixed_file=None, synth=False, includes=None):
  """
  Start compiling a specification.
  
  Currently, it produces a .des file that can be used by sequence designers.
  """
  
  print "Compiling '%s' ..." % basename
  # Read in system (or component)
  system = load_file(basename, args, prefix="", includes=includes)
  
  if fixed_file:
    print "Fixing sequences from file '%s'" % fixed_file
    fixed_sequences = load_fixed(fixed_file)
    for type_, name, fixed_seq in fixed_sequences:
      if type_ in "sequence":
        system.seqs[name].fix_seq( fixed_seq )
      elif type_ in "signal":
        # As a small hack, fix the first sequence in the list for the signal.
        for seq in system.signals[name]:
          if not seq[2]:
            seq[0].fix_seq( fixed_seq )
          else:
            seq[0].wc.fix_seq( fixed_seq )
      elif type_ == "strand":
        system.strands[name].fix_seq( fixed_seq )
      elif type_ == "structure":
        system.structs[name].fix_seq( fixed_seq )
  

  # Write the Zadeh-style design file
  print "System/component compiled into '%s'" % outputname
  outfile = open(outputname, "w")
  outfile.write("## Specification for %s compiled at: %s\n" % (basename, time.ctime()))
  if synth:
    system.output_synthesis("", outfile)
  else:
    system.output_nupack("", outfile)
  outfile.close()

  # Save compiler state to be reloaded when designer finishes
  print "Compiler state saved into '%s'" % savename
  save(system, savename)
  print "Run a designer on '%s' and process the result with python finish.py" % outputname

def save(obj, filename):
  """Save an object for later finishing."""
  f = open(filename, "w")
  pickle.dump(obj, f)
  f.close()

def load(filename):
  """Load it back."""
  if not os.path.isfile(filename):
    error("Cannot load save state. No such file '%s'." % filename)
  f = open(filename, "r")
  obj = pickle.load(f)
  f.close()
  return obj

if __name__ == "__main__":
  import sys
  from optparse import OptionParser
  
  try:
    import config_choices
  except ImportError:
    warning("DNA Circuit Compiler is not configured, please run config.py")
  
  # Parse command line options.
  usage = "usage: %prog [options] BASENAME [parameters ...]"
  parser = OptionParser(usage=usage)
  parser.set_defaults(pil=True) #verbose=True)
  #parser.add_option("-v", "--verbose", action="store_true", dest="verbose")
  # TODO: implement quiet
  #parser.add_option("-q", "--quiet", action="store_false", dest="verbose")
  parser.add_option("--fixed", help="Fix specific sequences listed in FILE", metavar="FILE")
  parser.add_option("--pil", action="store_true", help="Output in the new .pil format [Default]")
  parser.add_option("--des", action="store_false", dest="pil", help="Output in Zadeh's .des format instead of .pil format")
  parser.add_option("--synthesis", action="store_true", dest="pil", help="Depricated, use --pil instead.")
  parser.add_option("--output", help="Output file [defaults to BASENAME.pil]", metavar="FILE")
  parser.add_option("--save", help="Saved state file [defaults to BASENAME.save]", metavar="FILE")
  parser.add_option("-I","--include", help="Add PATH to search path for component imports", metavar="PATH",action="append")
  (options, args) = parser.parse_args()
  
  # Get basename of input specification
  if len(args) < 1:
    parser.error("missing required argument BASENAME")
  basename = args[0]
  # Infer the basename if a full filename is given
  p = re.match(r"(.*)\.(sys|comp)\Z", basename)
  if p:
    basename = p.group(1)
  
  # Set filename defaults
  if not options.output:
    if options.pil:
      options.output = basename + ".pil"
    else:
      options.output = basename + ".des"
  if not options.save:
    options.save = basename + ".save"
  
  # Eval remaining arguments
  import quickargs
  args, keys = quickargs.get_args(args[1:])
  assert not keys, "Don't provide keywords to compiler.py"
  
  compiler(basename, args, options.output, options.save, options.fixed, options.pil, options.include)
