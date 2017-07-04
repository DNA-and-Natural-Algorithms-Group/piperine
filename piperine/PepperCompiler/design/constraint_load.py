import string
from collections import defaultdict

from constraints import propagate_constraints
from PIL_parser import load_spec
from PIL_DNA_classes import group, rev_group, complement, seq_comp

from ..utils import warning

def min_(xs):
  """Returns min element (or None if there are no elements)"""
  if len(xs) == 0:
    return None
  else:
    return min(xs)

def intersect_groups(x1, x2):
  g1 = group[x1]; g2 = group[x2]
  inter = set(g1).intersection(set(g2))
  # TODO: better errors for overconstrained system
  if len(inter) == 0:
      raise ValueError
  inter = list(inter)
  inter.sort()
  inter = string.join(inter, "")
  return rev_group[inter]

class Constraints(object):
  def __init__(self):
    self.eq = {}
    self.wc = {}
    self.st = {}
    self.name = {}
  
  def init(self, x, letter="N", name="NoName"):
    """Must initialize an index 'x' before constraining it."""
    assert x not in self.eq and x not in self.wc and x not in self.st, "Cannot initialize an index twice: %r" % x
    self.eq[x] = []
    self.wc[x] = []
    self.st[x] = letter
    self.name[x] = name
  
  def add_eq(self, x, y):
    """Add a single equality constraint between items x and y."""
    self.eq[x].append(y)
    self.eq[y].append(x)
  
  def add_eqs(self, x, y, num):
    """Add equality constraints for 'num' bases starting with x and y"""
    for i in range(num):
      self.add_eq(x + i, y + i)
  
  def add_wc(self, x, y):
    """Add a single complementarity constraint between items x and y."""
    self.wc[x].append(y)
    self.wc[y].append(x)
  
  def add_wcs(self, x, y, num):
    """Add complementarity constraints for 'num' bases starting with x and y"""
    for i in range(num):
      self.add_wc(x + i, y - i)

  def propagate(self):
    self.eq, self.wc = propagate_constraints(self.eq, self.wc)
  
  def propagate_templates(self):
    """
    Propagate the constraints on sequence specified by the templates of equal or 
    complementary nucleotides.
    
    Warning: Only run after propagating the eq and wc constraints.
    """
    assert set(self.st.keys()) == set(self.eq.keys()) == set(self.wc.keys())
    done = set()
    for x, st in self.st.items():
      if x not in done:
        # Constraint must match all equal ...
        for y in self.eq[x]:
          assert y not in done, (x, y)
          st_y = self.st[y]
          try:
              st = intersect_groups(st, st_y)
          except ValueError:
              raise ValueError("{}, {}".format( self.name[x], self.name[y] ) )
        # ... and the complement of all wc
        for y in self.wc[x]:
          assert y not in done, (x, y)
          st_y = complement[self.st[y]]
          try:
              st = intersect_groups(st, st_y)
          except ValueError:
              raise ValueError("{}, {}".format( self.name[x], self.name[y] ) )
        
        # Apply new template
        for y in self.eq[x]:
          self.st[y] = st
          done.add(y)
        for y in self.wc[x]:
          self.st[y] = complement[st]
          done.add(y)
  
  def get_reps(self, isvalid=(lambda y: isinstance(y, int))):
    """Get representatives for eq and wc classes. Reps are minimum valid elements"""
    assert set(self.eq.keys()) == set(self.wc.keys())
    self.eq_rep = {}
    self.wc_rep = {}
    
    for x in self.eq.keys():
      self.eq_rep[x] = min_([y for y in self.eq[x] if isvalid(y)])
      self.wc_rep[x] = min_([y for y in self.wc[x] if isvalid(y)])
      for y in self.eq[x]:
        self.eq_rep[y] = self.eq_rep[x]
        self.wc_rep[y] = self.wc_rep[x]
      for y in self.wc[x]:
        self.eq_rep[y] = self.wc_rep[x]
        self.wc_rep[y] = self.eq_rep[x]
  
  def dump(self):
    """
    Output st constraints and eq and wc representatives as standard lists.
    These are not quite in spuriousC form because they still
      use None and base 0 indexing instead of 
      using -1 and base 1 indexing
    """
    assert set(self.st.keys()) == set(self.eq.keys()) == set(self.wc.keys())
    ## NOTE: drops the sequences and super-sequences from constraints lists.
    N = max([key for key in self.eq_rep.keys() if isinstance(key, int)]) + 1 # Number of indexes to be used
    eq = [self.eq_rep.get(i) for i in range(N)]
    wc = [self.wc_rep.get(i) for i in range(N)]
    st = [self.st.get(i) for i in range(N)]
    
    return eq, wc, st


def index_func_struct(spec):
  """Return an index function based on the STRUCTURE-centric specification"""
  struct_start = [None] * len(spec.structs)
  strand_start = [None] * len(spec.strands)
  for num, strand in enumerate(spec.strands.values()):
    strand.num = num
  
  prev_length = 0  # Keep a running sum of the lengths of structures
  for num, struct in enumerate(spec.structs.values()):
    struct.num = num
    struct_start[num] = prev_length
    for strand in struct.strands:
      if strand_start[strand.num] == None:
        strand_start[strand.num] = prev_length
      # On the linear constraint array we will put one blanks between strands
      prev_length += strand.length + 1
    # and two blanks between structures, 
    prev_length += 1
  
  def get_index(struct, index):
    """Return the general index for a nucleotide from the given structure."""
    assert index < struct.length, "get_index called with an index that is too big. %d >= %d" % (index, struct.length)
    result = struct_start[struct.num]
    for strand in struct.strands:
      # If this index puts us past this strand
      if index >= strand.length:
        index -= strand.length
        result += strand.length + 1
      # Otherwise it leaves us midway in this strand
      else:
        return result + index
    # We should never finish this loop (only happen if struct.length > sum(strands lengths)
    assert False, "Structure %s has length %d, but strands total length %d" % (struct.name, struct.length, sum([strand.length for strand in struct.strands]))
  
  def get_index_strand(strand, index):
    """Return the general index for a nucleotide from the given strand."""
    assert index < strand.length, "get_index called with an index that is too big. %d >= %d" % (index, strand.length)
    return strand_start[strand.num] + index

  def get_struct(index):
    """Return the structure for a given index."""
    raise Exception("Not Yet Implemented")

  def get_strand(index):
    """Return the strand and offset for a given index."""
    raise Exception("Not Yet Implemented")

  return get_index, get_index_strand, get_struct, get_strand

def index_func_strand(spec):
  """Return an index function based on the STRAND-centric specification"""
  ### TODO: Deal with dummy strands
  strand_start = [None] * len(spec.strands)
  prev_length = 0  # Keep a running sum of the lengths of strands
  for num, strand in enumerate(spec.strands.values()):
    strand.num = num
    strand_start[num] = prev_length
    # On the linear constraint array we will put two blanks between strands
    prev_length += strand.length + 2
  
  def get_index_strand(strand, index):
    """Return the general index for a nucleotide from the given strand."""
    assert index < strand.length, "get_index called with an index that is too big. %d >= %d" % (index, strand.length)
    return strand_start[strand.num] + index
  
  def get_index(struct, index):
    """Return the general index for a nucleotide from the given structure."""
    assert index < struct.length, "get_index called with an index that is too big. %d >= %d" % (index, struct.length)
    for strand in struct.strands:
      # If this index puts us past this strand
      if index >= strand.length:
        index -= strand.length
      # Otherwise it leaves us midway in this strand
      else:
        return get_index_strand(strand, index)
    # We should never finish this loop (only happen if struct.length > sum(strands lengths)
    assert False, "Structure %s has length %d, but strands total length %d" % (struct.name, struct.length, sum([strand.length for strand in struct.strands]))

  def get_struct(index):
    """Return the structure for a given index."""
    raise Exception("Not Yet Implemented")

  def get_strand(index):
    """Return the strand and offset for a given index."""
    if type(index) == tuple:
        index = index[0]
    num = len([x for x in strand_start if index > x])-1
    return (spec.strands.values()[num], index - strand_start[num])

  return get_index, get_index_strand, get_struct, get_strand

class Convert(object):
  def __init__(self, filename, struct_orient=False):
    """
    Load a file (in PIL format) and prepare to get constraints.
    
    struct_orient == True if we are listing out all structures in constraint files
    struct_orient == False if we are listing out strands (only once each)
    """
    ## Load the specification
    self.spec = load_spec(filename)
    
    ## Create constraint object
    self.constraints = Constraints()
    
    self.struct_orient = struct_orient
  
  def get_constraints(self):
    """Distil out constraints information from the specification."""
    ## Initialize all struct/strand constraints lists
    if self.struct_orient:
      self.get_index, self.get_index_strand, self.get_struct, self.get_strand = index_func_struct(self.spec)
      
      for struct in self.spec.structs.values():
        for x in range(struct.length):
          x2 = self.get_index(struct, x)
          self.constraints.init(x2)
      
      # Constrain all instances of the same strand to be equal
      for struct in self.spec.structs.values():
        offset = 0
        for strand in struct.strands:
          for x in range(strand.length):
            x2 = self.get_index_strand(strand, x)
            y2 = self.get_index(struct, offset + x)
            self.constraints.add_eq(x2, y2)
          offset += strand.length
    
    else: # if strand oriented
      self.get_index, self.get_index_strand, self.get_struct, self.get_strand = index_func_strand(self.spec)
      
      for strand in self.spec.strands.values():
        for x in range(strand.length):
          x2 = self.get_index_strand(strand, x)
          self.constraints.init(x2, name=strand.name)
    
    self.constraints.get_strand = self.get_strand 
    ## Add constraints
    # Structural constraints
    for struct in self.spec.structs.values():
      for x, y in struct.bonds:
        x2 = self.get_index(struct, x)
        y2 = self.get_index(struct, y)
        self.constraints.add_wc(x2, y2)
    
    
    ## Initialize all sequence constraints lists
    num = 0
    for seq in self.spec.base_seqs.values():
      seq.num = num
      num += 1
      for x, letter in enumerate(seq.template):
        x2 = (seq.num, x)
        self.constraints.init(x2, letter, name=seq.name)
      # We list both sequences and reverse_sequences making constraints code simpler
      seq.wc.num = num
      num += 1
      for x in range(seq.wc.length):
        x2 = (seq.wc.num, x)
        self.constraints.init(x2, name=seq.name+"_reversed")
        self.constraints.add_wc(x2, (seq.num, seq.length - x - 1) ) # Add wc constraint
    # and for super-sequences
    for seq in self.spec.sup_seqs.values():
      seq.num = num
      num += 1
      for x in range(seq.length):
        x2 = (seq.num, x)
        self.constraints.init(x2)
      # and wc
      seq.wc.num = num
      num += 1
      for x in range(seq.wc.length):
        x2 = (seq.wc.num, x)
        self.constraints.init(x2, name=seq.name+"_super")
        self.constraints.add_wc(x2, (seq.num, seq.length - x - 1) ) # Add wc constraint
    
    
    # Equality constraints
    for eqlist in self.spec.equals:
      first_seq = eqlist[0]
      for seq in eqlist[1:]:
        assert seq.length == first_seq.length
        for x in range(seq.length):
          self.constraints.add_eq( (first_seq.num, x), (seq.num, x) )
    
    # Super-sequence constraints
    for seq in self.spec.sup_seqs.values():
      offset = 0
      for sub_seq in seq.seqs:
        for x in range(sub_seq.length):
          self.constraints.add_eq( (seq.num, offset + x), (sub_seq.num, x) )
        offset += sub_seq.length
    
    # Strand constraints
    for strand in self.spec.strands.values():
      offset = 0
      for seq in strand.seqs:
        for x in range(seq.length):
          x2 = self.get_index_strand(strand, offset + x)
          self.constraints.add_eq(x2, (seq.num, x) )
        offset += seq.length
    
    ## Propagate constraints
    self.constraints.propagate()
    
    ## Propagate templating constraints
    self.constraints.propagate_templates()
    
    ## Get constraints
    self.constraints.get_reps()
    eq, wc, st = self.constraints.dump()
    
    return eq, wc, st

  def process_results(self, nts):
    """Once a designer has designed a nucleotide sequence, reincorporate that info back into the specification."""
    for strand in self.spec.strands.values():
      seq = [nts[self.get_index_strand(strand, x)] for x in range(strand.length)]
      strand.set_seq(string.join(seq, ""))
    
    for struct in self.spec.structs.values():
      struct.get_seq()
  
  def output(self, outname, findmfe=True):
    """Output designed sequences to outfilename in .mfe format."""
    f = open(outname, "w")
    for num, struct in enumerate(self.spec.structs.values()):
      # TODO-maybe: change output format so we don't need to run DNAfold in designer.
      if findmfe:
        from ..DNAfold import DNAfold
        struct.mfe_struct, dG = DNAfold(struct.seq)
      else:
        struct.mfe_struct = struct.struct
        dG = 0
      
      # Write structure (with dummy content)
      f.write("%d:%s\n" % (num, struct.name))
      gc_content = (struct.seq.count("C") + struct.seq.count("G")) / struct.length
      f.write("%s %f %f %d\n" % (struct.seq, 0, gc_content, 0))
      f.write("%s\n" % struct.struct)   # Target structure
      f.write("%s\n" % struct.mfe_struct)  # MFE structure
    
    for seq_num, seq in enumerate(self.spec.seqs.values()):
      num = seq_num + len(self.spec.structs)
      
      # Deal with sequences that haven't been designed (because they were not used in strands).
      if not seq.seq:
        seq.seq = seq.get_seq()
        seq.wc.seq = seq_comp(seq.seq)
        # TODO: Find some way to deal with this appropriately.
        # Right now we just return NNNN...
        if "N" in seq.seq:
          pass #warning("Sequence %s was not designed because it was never used in a strand." % seq.name)
      
      # Write sequence (with dummy content)
      f.write("%d:%s\n" % (num, seq.name))
      gc_content = (seq.seq.count("C") + seq.seq.count("G")) / seq.length
      f.write("%s %f %f %d\n" % (seq.seq, 0, gc_content, 0))
      f.write(("."*seq.length+"\n")*2) # Write out dummy structures.
      # Write wc of sequence (with dummy content)
      seq = seq.wc
      f.write("%d:%s\n" % (0, seq.name))
      f.write("%s %f %f %d\n" % (seq.seq, 0, gc_content, 0))
      f.write(("."*seq.length+"\n")*2) # Write out dummy structures.
    
    f.write("Total n(s*) = %f" % 0) # Dummy n(s*)
    f.close()

if __name__ == "__main__":
  import sys
  
  convert = Convert(sys.argv[1], (len(sys.argv) > 2))
  print convert.get_constraints()
