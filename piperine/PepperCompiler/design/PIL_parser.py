import re
import sys

import PIL_class
from DNA_nupack_classes import group

from ..utils import error, match

def load_spec(filename):
  """Load a PIL style input specification."""
  f = open(filename, "r")
  
  # Create new specification object ...
  spec = PIL_class.Spec()
  # ... and populate it with the contents of the file.
  for line in f:
    # Strip comments and whitespace
    line = re.sub(r"#.*\n", "", line)
    line = line.strip()
    # Skip empty lines
    if not line:
      continue
    
    #print line,
    # Read the command name off (if there is a command name)
    command = line.split()[0]
    if command == "sequence":
      name, template = parse_seq(line)
      spec.add_seq(name, template)
    elif command in ("super-sequence", "sup-sequence"):
      name, sub_seq_names = parse_sup_seq(line)
      spec.add_sup_seq(name, sub_seq_names)
    elif command == "strand":
      name, seq_names, dummy = parse_strand(line)
      spec.add_strand(name, seq_names, dummy)
    elif command == "structure":
      name, strand_names, struct, params = parse_struct(line)
      spec.add_struct(name, strand_names, struct, params)
    elif command == "equal":
      seq_names = parse_equal(line)
      spec.add_equal(seq_names)
    elif command == "kinetic":
      pass
    else:
      error("Parse Error in file '%s': Command '%s' not valid.\n%s" % (filename, command, line)) # TODO: more info
  
  return spec

def parse_seq(line):
  """Parse sequence statements"""
  m = match(r"sequence ([\w-]+) = ([^:\s]*)( : .*)?", line)
  if not m:
    error("Invalid sequence statement format:\n"
          "Should be: sequence <name> = <constraint template>\n"
          "Was:       %s" % line)
  name, template = m.group(1, 2)
  if not set(template).issubset( set(group.keys()) ):
    error("Sequence's constraint template must be in allowed alphabet (%r).\nLine: %s" % (group.keys(), line))
  return name, template

def parse_sup_seq(line):
  """Parse super-sequence statements"""
  m = match(r"sup(er)?-sequence ([\w-]+) = ([^:]*)( : .*)?", line)
  if not m:
    error("Invalid super-sequence statement format:\n"
          "Should be: super-sequence <name> = <sub-sequence names>\n"
          "Was:       %s" % line)
  name, seq_names = m.group(2, 3)
  seq_names = seq_names.split()
  return name, seq_names

def parse_strand(line):
  """Parse strand statements"""
  m = match(r"strand (\[dummy\] )?([\w-]+) = ([^:]*)( : .*)?", line)
  if not m:
    error("Invalid strand statement format:\n"
          "Should be: strand <name> = <sequence names>\n"
          "or:        strand [dummy] <name> = <sequence names>\n"
          "Was:       %s" % line)
  dummy, name, seq_names = m.group(1, 2, 3)
  dummy = bool(dummy)
  seq_names = seq_names.split()
  return name, seq_names, dummy

def parse_struct(line):
  """Parse structure statements"""
  m = match(r"structure (\[(\w+)\])? ([\w-]+) = ([^:]*) : (.*)", line)
  if not m:
    error("Invalid structure statement format:\n"
          "Should be: structure <name> = <strand names> : <secondary structure>\n"
          "or:        structure [<parameters>] <name> = <strand names> : <secondary structure>\n"
          "Was:       %s" % line)
  params, name, strand_names, struct = m.group(2, 3, 4, 5)
  strand_names = strand_names.split("+")
  strand_names = [sname.strip() for sname in strand_names]  # Clean off whitespace
  struct = struct.replace(" ", "").replace("\t", "") # Clear out whitespace
  if not set(struct).issubset( set(".()+") ):
    error('Secondary structure must only use allowd alphabet (".()+").\nLine: %s' % line)
  return name, strand_names, struct, params

def parse_equal(line):
  """Parse super-sequence statements"""
  seq_names = line.split()[1:]
  return seq_names

if __name__ == "__main__":
  import sys
  
  print load_spec(sys.argv[1])
