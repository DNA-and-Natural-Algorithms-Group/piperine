import re
import sys


from ..utils import error
import nupack_in_class
from DNA_nupack_classes import group

def load_file(filename):
  """Load a NUPACK style input file."""
  f = open(filename, "r")
  
  # Create new specification object ...
  spec = nupack_in_class.Spec()
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
    pieces = line.split() # TODO: Put in try/except for when this fails
    if pieces[0] == "structure":
      name, struct = parse_struct(line)
      spec.add_structure(name, struct)
    elif pieces[0] == "sequence":
      name, seq = parse_seq(line)
      spec.add_sequence(name, seq)
    elif pieces[0] == "prevent":
      pass
    elif pieces[1] == ":": # Applying sequences to a structure
      struct_name, seqs = parse_apply(line)
      spec.add_apply(struct_name, seqs)
    elif pieces[1] == "<": # Set the objective function for a structure
      pass
    else:
      error("Parse Error in file '%s': Command '%s' not valid.\n%s" % (filename, pieces[0], line)) # TODO: more info
  
  return spec

def parse_struct(line):
  command, name, eq, struct = line.split(None, 3) # TODO: Put in try/except for when this fails
  assert eq == "=", "Structure syntax incorrect" # TODO: add line and syntax
  return name, struct

def parse_seq(line):
  try:
    sequence, name, eq, seq = line.split(None, 3) # TODO: Put in try/except for when this fails
    assert eq == "=", "Sequence syntax incorrect" # TODO: add line and syntax
  except (ValueError, AssertionError):
    error("Invalid sequence statement format:\n"
          "Should be: sequence <name> = <constraints>\n"
          "Invalid:   %s" % line)
  assert set(seq).issubset( set(group.keys()) ), "Sequence constraints must be written out in allowed alphabet. %r not in %r" % (seq, group.keys())
  return name, seq

def parse_apply(line):
  struct_name, colon, seqs = line.split(None, 2)
  
  temp = seqs.split()
  seqs = []
  for item in temp:
    p = re.match(r"([\w_-]+)(\*)?", item)
    seq_name, wc = p.groups("")     # wc = "*" if this is the complement, or "" otherwise
    seqs.append((seq_name, wc))
  
  return struct_name, seqs

if __name__ == "__main__":
  import sys
  
  load_file(sys.argv[1])
