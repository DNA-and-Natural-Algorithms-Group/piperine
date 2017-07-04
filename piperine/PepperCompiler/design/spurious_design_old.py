#!/usr/bin/env python
"""
Designs sequences using Winfree's SpuriousDesign/spuriousC.c algorithm.
Uses Joe Zadah's input and output formats for compatibility with compiler.
"""

import string
import subprocess
import os

from new_loading import load_file
from DNA_nupack_classes import group, rev_group, complement, seq_comp

from ..DNAfold import DNAfold
from ..utils import error

DEBUG = False

# HACK
group["_"] = ""
rev_group[""] = "_"
complement["_"] = "_"

def intersect_groups(x1, x2):
  g1 = group[x1]; g2 = group[x2]
  inter = set(g1).intersection(set(g2))
  #assert len(inter) > 0, "System over-constrained. %s and %s cannot be equal." % (x1, x2)
  inter = list(inter)
  inter.sort()
  inter = string.join(inter, "")
  return rev_group[inter]

# SpuriousC notation for nothing (specifically, there being no wc complement).
NOTHING = -1

def min_(foo):
  """Return the min element (or NOTHING if there are no elements)."""
  if len(foo) == 0:
    return NOTHING
  else:
    return min(foo)

def min_rep(*terms):
  """Return the minimum representative. Where NOTHING is no representative."""
  terms = [x for x in terms if x != NOTHING]
  return min_(terms)

class Connections(object):
  """
  Keeps track of which parts of dna sequence are constrained to be equal or 
  complementary to which other parts.
  """
  def __init__(self, spec):
    # TODO: split strands in a structure.
    
    # seq_const[n][0] = list of (s, i) where spec.structs[s].seq[m] is spec.seqs[n]
    # seq_const[n][1] = list of (s, i) where it's the complement
    seq_const = [([], []) for seq in spec.seqs]
    
    self.structs = spec.structs.values()
    self.seqs = spec.seqs.values()
    
    # HACK: Adding pseudo-structures for sequences
    self.num_structs = len(self.structs) + len(self.seqs)
    
    self.struct_length = [None] * self.num_structs
    
    # Create seq_const and find lengths of each structure.
    for s, struct in enumerate(self.structs):
      i = 0
      for seq in struct.seqs:
        # Add this struct, index to the appropriate set.
        seq_const[seq.num][seq.reversed].append((s, i))
        i += seq.length
      self.struct_length[s] = i
    # And for each sequence pseudo-structure
    for seq_num, seq in enumerate(self.seqs):
      s = len(self.structs) + seq_num
      seq_const[seq.num][seq.reversed].append((s, 0))
      self.struct_length[s] = seq.length
    
    self.table = {NOTHING: (NOTHING, NOTHING)}
    # Expand these sets for specific nucleotides ...
    for seq_num, (eq, wc) in enumerate(seq_const):
      seq = spec.seqs.get_index(seq_num)
      for nt in xrange(seq.length):
        new_eq = [(s, i + nt) for (s, i) in eq]
        rep_eq = min_(new_eq)
        new_wc = [(s, i + (seq.length - 1 - nt)) for (s, i) in wc]
        rep_wc = min_(new_wc)
        # ... and save those into the connection table.
        for (s, i) in new_eq:
          self.table[(s, i)] = (rep_eq, rep_wc)
        for (s, i) in new_wc:
          self.table[(s, i)] = (rep_wc, rep_eq)
  
  def apply_wc(self, x, y):
    """Apply WC complementarity condition between items."""
    eq_x, wc_x = self.table[x]
    eq_y, wc_y = self.table[y]
    # x is complementary to y so merge their representatives.
    new_eq_x = new_wc_y = min_rep(eq_x, wc_y)
    new_eq_y = new_wc_x = min_rep(eq_y, wc_x)
    
    for z, (eq_z, wc_z) in self.table.items():
      if eq_z != NOTHING:
        # If z == x (== y*), update it to be equal
        if eq_z in (eq_x, wc_y):
          self.table[z] = new_eq_x, new_wc_x
        # If z == y (== x*), update it
        elif eq_z in (eq_y, wc_x):
          self.table[z] = new_eq_y, new_wc_y
  
  def printf(self):
    print "Connections.printf()"
    for s in xrange(self.num_structs):
      for i in xrange(self.struct_length[s]):
        (eq, wc) = self.table[(s, i)]
        print (s, i), eq, wc
    print

def prepare(in_name):
  """Create eq, wc and sequence template lists in spuriousC format."""
  # Load specification from input file.
  spec = load_file(in_name)
  
  # Load structures and subsequences into Connections object
  # and connect subsequences that are equal and complementary.
  c = Connections(spec)
  
  # Connect regions that are constrained to helixes.
  for s, struct in enumerate(c.structs):
    for start, stop in struct.bonds:
      c.apply_wc((s, start), (s, stop))
  
  # Convert from (strand, index) format to multi_index format
  # Create the conversion function
  f = {NOTHING: NOTHING}
  eq = []
  wc = []
  st = []  # We start it as a list because python strings aren't mutable
  for s, struct in enumerate(c.structs):
    start = len(eq)
    for i in xrange(c.struct_length[s]):
      this = len(eq)
      if struct.struct[this-start] == "+":
        # Single spaces between strands.
        eq += [NOTHING]
        wc += [NOTHING]
        st += [" "]
      f[(s, i)] = len(eq)
      # Get the representatives for (s, i)
      eq_si, wc_si = c.table[(s, i)]
      # and convert them
      eq.append(eq_si)
      wc.append(wc_si)
      st.append(struct.seq[i])
    # Double spaces between structures
    eq += [NOTHING, NOTHING]
    wc += [NOTHING, NOTHING]
    st += [" ", " "]
  # Mark the end of the real structures (next are pseudo-structures we will remove later
  real_length = len(st) - 2
  assert len(st) == len(eq) == len(wc)
  # The pseudo-structures.
  for seq_num, seq in enumerate(c.seqs):
    s = len(c.structs) + seq_num
    for i in xrange(c.struct_length[s]):
      f[(s, i)] = len(eq)
      # Get the representatives for (s, i)
      eq_si, wc_si = c.table[(s, i)]
      # and convert them
      eq.append(eq_si)
      wc.append(wc_si)
    # Double spaces between structures
    eq += [NOTHING, NOTHING]
    wc += [NOTHING, NOTHING]
    st += seq.seq
    st += [" ", " "]
    
  
  # Update using the conversion function
  eq = [f[x] for x in eq[:-2]]  # eq[:-2] to get rid of the [NOTHING, NOTHING] at the end
  wc = [f[x] for x in wc[:-2]]
  st = st[:-2]
  
  # Constrain st appropriately
  for i in xrange(len(eq)):
    #sys.stdout.write(st[i])
    # Skip strand breaks
    if eq[i] == NOTHING:
      continue
    if eq[i] < i:
      j = eq[i]
      old_stj = st[j]
      st[j] = intersect_groups(st[j], st[i])
      if st[j] == "_" and "_" not in (st[i], old_stj):
        print
        print i, j, st[i], old_stj
        
    if wc[i] != NOTHING and wc[i] < i:
      j = wc[i]
      old_stj = st[j]
      st[j] = intersect_groups(st[j], complement[st[i]])
      if st[j] == "_" and "_" not in (st[i], old_stj):
        print
        print i, j, st[i], old_stj
  # Propagate the changes
  for i in xrange(len(eq)):
    if eq[i] != NOTHING and eq[i] < i:
      j = eq[i]
      st[i] = st[j]
    if wc[i] != NOTHING:
      j = wc[i]
      st[i] = complement[st[j]]
  
  # Get rid of pseudo structures
  st = st[:real_length]
  # Also clear any pointers to these pseudo-structures
  eq = [x if x < real_length else NOTHING for x in eq[:real_length]]
  wc = [x if x < real_length else NOTHING for x in wc[:real_length]]
  return st, eq, wc, c

def process_result(c, inname, outname):
  """Output sequences in Joe's format."""
  # Read spuriousC's output.
  f = open(inname, "r")
  nts = f.read()
  f.close()
  # File has all sorts of runtime info.
  # The final sequences are stored on the last full line.
  nts = nts.split("\n")[-2]
  #print repr(nts)
  seqs = string.split(nts, "  ") # Sequences for each structure
  
  f = open(outname, "w")
  assert len(c.structs) == len(seqs), "%d != %d" % (len(c.structs), len(seqs))
  for struct, nseq in zip(c.structs, seqs):
    if DEBUG: print struct.name, nseq
    # Save the sequence
    struct.nseq = nseq.replace(" ", "+")
    struct.mfe_struct, dG = DNAfold(struct.nseq)
    #print repr(struct.nseq)
    
    # Extract sequences for domains (seqs)
    full_seq = struct.nseq.replace("+", "")
    i = 0 # index of start of next domain
    for seq in struct.seqs:
      seq.nseq = full_seq[i:i+seq.length]
      seq.wc.nseq = seq_comp(seq.nseq)
      i += seq.length
    assert i == struct.length, "Sequences in struct %s have total length %d != %d, the length of the struct" % (struct.name, i, struct.length)
    
    # Write structure (with dummy content)
    f.write("%d:%s\n" % (0, struct.name))
    gc_content = (struct.nseq.count("C") + struct.nseq.count("G")) / struct.length
    f.write("%s %f %f %d\n" % (struct.nseq, 0, gc_content, 0))
    f.write("%s\n" % struct.struct)   # Target structure
    f.write("%s\n" % struct.mfe_struct)  # Dummy MFE structure
  
  for seq_num, seq in enumerate(c.seqs):
    # TODO: not returning some sequences may crash finish.py
    if seq.nseq:
      # Write sequence (with dummy content)
      f.write("%d:%s\n" % (0, seq.name))
      gc_content = (seq.nseq.count("C") + seq.nseq.count("G")) / seq.length
      f.write("%s %f %f %d\n" % (seq.nseq, 0, gc_content, 0))
      f.write(("."*seq.length+"\n")*2) # Write out dummy structures.
      # Write wc of sequence (with dummy content)
      seq = seq.wc
      f.write("%d:%s\n" % (0, seq.name))
      #wc_seq = string.join([complement[symb] for symb in seq.nseq[::-1]], "")
      f.write("%s %f %f %d\n" % (seq.nseq, 0, gc_content, 0))
      f.write(("."*seq.length+"\n")*2) # Write out dummy structures.
  f.write("Total n(s*) = %f" % 0)
  f.close()

def print_list(xs, filename, format):
  """Prints a list 'xs' to a file using space separation format."""
  f = open(filename, "w")
  for x in xs:
    f.write(format % x)
  f.close()

def design(basename, infilename, outfilename, cleanup, verbose=False, reuse=False, extra_pars=""):
  stname = basename + ".st"
  wcname = basename + ".wc"
  eqname = basename + ".eq"
  sp_outname = basename + ".sp"
  
  #if reuse:
  #  for name in stname, wcname, eqname:
  #    assert os.path.isfile(name), "Error: requested --reuse, but file '%s' doesn't exist" % name
  
  # Prepare the constraints
  print "Reading design from  file '%s'" % infilename
  print "Preparing constraints files for spuriousC."
  st, eq, wc, c = prepare(infilename)
  
  # Convert specifications
  eq = [x+1 for x in eq]
  for i, x in enumerate(wc):
    if x != -1:
      wc[i] += 1
  
  # Print them to files
  print_list(st, stname, "%c")
  print_list(wc, wcname, "%d ")
  print_list(eq, eqname, "%d ")
  
  if "_" in st:
    print "System over-constrained."
    sys.exit(1)
  
  # Run SpuriousC
  # TODO: take care of prevents.
  if verbose:
    quiet = "quiet=ALL | tee %s" % sp_outname
  else:
    quiet = "quiet=TRUE > %s" % sp_outname
  
  command = "spuriousC score=automatic template=%s wc=%s eq=%s %s %s" % (stname, wcname, eqname, extra_pars, quiet)
  print command
  subprocess.check_call(command, shell=True)
  
  # Process results
  print "Processing results of spuriousC."
  process_result(c, sp_outname, outfilename)
  print "Done, results saved to '%s'" % outfilename
  if cleanup:
    print "Deleting temporary files"
    os.remove(stname)
    os.remove(wcname)
    os.remove(eqname)
    os.remove(sp_outname)

if __name__ == "__main__":
  import re
  from optparse import OptionParser
  
  from find_file import find_file, BadFilename
  
  if sys.version_info < (2, 5):
    error("Must use python 2.5 or greater.")
    
  # Parse command line options.
  usage = "usage: %prog [options] infilename [spuriousC_parameters ...]"
  parser = OptionParser(usage=usage)
  parser.set_defaults(verbose=False, cleanup=True)
  parser.add_option("-v", "--verbose", action="store_true", dest="verbose", help="Verbose output from spuriousC")
  parser.add_option("-q", "--quiet", action="store_false", dest="verbose", help="No output from spuriousC [Default]")
  parser.add_option("-o", "--output", help="Output file [defaults to BASENAME.mfe]", metavar="FILE")
  
  parser.add_option("--keep-temp", action="store_false", dest="cleanup", help="Keep temporary files (.st, .wc, .eq, .sp)")
  parser.add_option("--cleanup", action="store_true", dest="cleanup", help="Remove temporary files after use [Default]")
  # TODO: parser.add_option("--reuse", action="store_true", help="Reuse the .st, .wc and .eq files if they already exist (Saves time if a session was terminated, or if you want to rerun a design).")
  (options, args) = parser.parse_args()
  options.reuse = False # TODO
  
  if len(args) < 1:
    parser.error("missing required argument infilename")
  try:
    infilename = find_file(args[0])
  except BadFilename:
    parser.error("File not found: neither %s nor %s.des exist. Please supply correct infilename." % (args[0], args[0]))
  
  # Infer the basename if a full filename is given
  basename = infilename
  p = re.match(r"(.*)\.des\Z", basename)
  if p:
    basename = p.group(1)
  
  # Set filename defaults
  if not options.output:
    options.output = basename + ".mfe"
  
  # Collect extra arguments for spuriousC
  spurious_pars = string.join(args[1:], " ")
  
  design(basename, infilename, options.output, options.cleanup, options.verbose, options.reuse, spurious_pars)
