"""DNA design container classes"""
import string

# Global DNA nt groups
group = {"A": "A", "T": "T", "C": "C", "G": "G",
         "W": "AT", "S": "CG", "M": "AC", "K": "GT", 
         "B": "CGT", "V": "ACG", "D": "AGT", "H": "ACT",
         "N": "ACGT"} # SpuriousC group codes
rev_group = dict([(v, k) for (k, v) in group.items()])  # A reverse lookup for group.
complement = {"A": "T", "T": "A", "C": "G", "G": "C",
              "W": "W", "S": "S", "M": "K", "K": "M",
              "B": "V", "V": "B", "D": "H", "H": "D",
              "N": "N"} # Should satisfy set(group[complement[X]]) == set(seq_comp(group[X]))
def seq_comp(seq):
  """Returns the WC complement of a nucleotide sequence."""
  return string.join([complement[nt] for nt in reversed(seq)], "")

class Sequence(object):
  """Container for sequences"""
  def __init__(self, name, template):
    self.name = name
    self.template = template  # The template of constraints on sequence
    self.seq = None # Sequence has not been designed yet
    self.length = len(self.template)
    self.reversed = False
    # Build the dummy sequence for the W-C complement
    self.wc = ReverseSequence(self)
  
  def set_seq(self, seq):
    """Set the sequence."""
    if self.seq: # If this sequence is already defined
      assert self.seq == seq, "Sequence %s was designed with 2 different sequences: %s and %s" % (self.name, self.seq, seq)
    else: # If it's not defined yet
      self.seq = seq
      self.wc.seq = seq_comp(seq)
  
  def get_seq(self):
    """Return designed sequence or template as default."""
    if self.seq:
      return self.seq
    else:
      return self.template
  
  def __repr__(self):
    return "Sequence(%(name)r, %(template)r)" % self.__dict__

class ReverseSequence(Sequence):
  """Complements of defined sequences"""
  def __init__(self, wc):
    self.name = wc.name + "*"
    self.template = seq_comp(wc.template)
    self.seq = None
    self.length = wc.length
    self.reversed = True
    self.wc = wc


class SuperSequence(object):
  """Logical grouping of sequences"""
  def __init__(self, name, sub_seqs):
    self.name = name
    self.seq = None # Stores the sequence once it has been designed
    self.seqs = sub_seqs
    
    # Find length and base_seqs
    self.base_seqs = []  # The atomic Sequence's that constrain this SuperSequence
    self.length = 0
    for sub_seq in self.seqs:
      if isinstance(sub_seq, Sequence):
        self.base_seqs.append(sub_seq)
      else:
        assert isinstance(sub_seq, SuperSequence)
        self.base_seqs += sub_seq.base_seqs
      self.length += sub_seq.length
    
    self.reversed = False
    self.wc = ReverseSuperSequence(self)
  
  def set_seq(self, seq):
    """Set the sequence and apply it to all subsequences."""
    assert len(seq) == self.length, "Designed sequence length mismatch. %d != %d" % (len(seq), self.length)
    if self.seq: # If this sequence is already defined
      assert self.seq == seq, "Sequence %s was designed with 2 different sequences: %s and %s" % (self.name, self.seqs, seqs)
    else: # If it's not defined yet
      self.seq = seq
      self.wc.seq = seq_comp(seq)
      i = 0
      for sub_seq in self.seqs:
        sub_seq.set_seq(seq[i:i+sub_seq.length])
        i += sub_seq.length
  
  def get_seq(self):
    """Return the sequence (or reconstruct it from subsequences)"""
    if self.seq:
      return self.seq
    else:
      return string.join([sub_seq.get_seq() for sub_seq in self.seqs], "")
  
  def __repr__(self):
    return "SuperSequence(%(name)r, %(seqs)r)" % self.__dict__

class ReverseSuperSequence(SuperSequence):
  def __init__(self, wc):
    self.name = wc.name + "*"
    self.seq = None # Stores the sequence once it has been designed
    self.seqs = [seq.wc for seq in reversed(wc.seqs)]
    self.base_seqs = [seq.wc for seq in reversed(wc.base_seqs)]
    self.length = wc.length
    self.reversed = True
    self.wc = wc


class Strand(SuperSequence):
  """Container for strands. Inherits from SuperSequence for convinience."""
  def __init__(self, name, seqs, dummy):
    SuperSequence.__init__(self, name, seqs)
    self.dummy = dummy
  def __repr__(self):
    return "Strand(%(name)r, %(seqs)r, %(dummy)r)" % self.__dict__


def get_bonds(struct):
  """Get a list of bonds in a dot-paren structure."""
  struct = struct.replace("+", "")  # Positions do NOT include strand breaks
  bonds = []  # Will be a list of pairs (open_postion, close_postion)
  open_pos = []  # A FILO of open parentheses positions
  
  for pos, symb in enumerate(struct):
    if symb == "(":
      open_pos.append(pos)
    elif symb == ")":
      assert len(open_pos) != 0
      start = open_pos.pop() # The most recent open-paren
      bonds.append( (start, pos) )
    else:
      assert symb == ".", "Structure '%s' not in dot-paren form" % struct
  
  return bonds

class Structure(object):
  """Container for structures/complexes"""
  def __init__(self, name, strands, struct, params):
    self.name = name
    self.params = params # Optimization parameters, TODO: actually deal with these
    self.struct = struct
    self.bonds = get_bonds(struct)
    self.strands = strands
    self.seq = None # Stores the sequence once it has been defined.
    
    # Find length and base_seqs
    self.base_seqs = []
    self.length = 0
    sub_structs = [strand_struct for strand_struct in self.struct.split("+")] # Check that lengths match up
    assert len(strands) == len(sub_structs), "Mismatch: Structure %s is defined by %d strands, but secondary structure has %d strands" % (name, len(strands), len(sub_structs))
    for strand, sub_struct in zip(strands, sub_structs):
      assert isinstance(strand, Strand), "Structure %s must get strands" % name
      assert strand.length == len(sub_struct), "Mismatch: Strand %s in structure %s has length %d, but sub-structure %s implies %d" % (strand.name, name, strand.length, sub_struct, len(sub_struct))
      self.base_seqs += strand.base_seqs
      self.length += strand.length
  
  def get_seq(self):
    """Get sequence from strands which have been set."""
    self.seq = string.join([strand.seq for strand in self.strands], "+")
  
  def __repr__(self):
    return "Structure(%(name)r, %(strands)r, %(struct)r, %(params)r)" % self.__dict__

