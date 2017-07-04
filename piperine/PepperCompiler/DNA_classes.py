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
              "N": "N"} # Should satisfy set(group[complement[X]]) == set(wc(group[X]))
def wc(seq):
  """Returns the WC complement of a nucleotide sequence."""
  return string.join([complement[nt] for nt in reversed(seq)], "")

class WildError(AssertionError):
  """For when a sequence is defined a wildcard, but without a length."""

WILDCARD = "?"

class Sequence(object):
  """Container for sequences"""
  def __init__(self, name, prefix, constraints, length=None):
    self.name = name
    self.prefix = prefix  # Just keep it for __repr__
    self.full_name = prefix + name
    self.seq = None # Stores the sequence once it has been defined.
    self.reversed = False
    
    self.length, self.const = self._get_length_const(length, list(constraints), name)
    
    # Sequences of length 0 must be treated specially, they are not passed on 
    #   to lower levels, but are allowed for the convinience of users.
    self.dummy = (self.length == 0)
    
    # Keep track of whether this sequence has been used in a strand yet.
    # If it never is, it will not be designed by most designers, so we warn user
    self.in_strand = False
    
    # Build the dummy sequence for the W-C complement
    self.wc = ReverseSequence(self)
  
  def _get_length_const(self, length, constraints, name):
    """Work out all the complicated wildcard rules to get the precise length and constraints list."""
    # Check length and resolve wildcard
    lengths = [num for num, code in constraints]
    wilds = lengths.count(WILDCARD)
    assert wilds <= 1, "Too many wildcards in sequence %s" % name
    if wilds == 0: # no wildcards
      seq_lengths = sum(lengths)
      if length:
        assert length == seq_lengths, "Length mismatch for sequence %s (%r != %r)" % (name, length, seq_lengths)
      else: # If length was not specified (None), we set it
        length = seq_lengths
    else: # one wildcard
      if length == None: raise WildError("Sequence %s has a ?. but no length specified" % name)
      check_length = sum([x for x in lengths if x != WILDCARD])
      wild_length = length - check_length  # Wildcard is set so that total length is right
      assert wild_length >= 0, "Sequence %s too short (%r > %r)" % (name, length, check_length)
      i = lengths.index(WILDCARD)
      constraints[i] = (wild_length, constraints[i][1])
    
    const = ""
    for (num, base) in constraints:
      const += base * num  # We represent constriants in long-form
    
    return length, const
  
  def fix_seq(self, fixed_seq):
    """Constrian ourselves to a specific sequence."""
    assert len(fixed_seq) == self.length
    for const_nt, fixed_nt in zip(self.const, fixed_seq):
      const_set = set(group[const_nt])
      fixed_set = set(group[fixed_nt])
      assert fixed_set.issubset( const_set ), "fix_seq: sequence %s is not a subset of constraint %s for sequence %s" % (fixed_seq, self.const, self.full_name)
    self.const = fixed_seq
  
  def __invert__(self):
    """Returns the Watson-Crick complementary sequence."""
    return self.wc
  def __repr__(self):
    return "Sequence(%(name)r, %(prefix)r, %(const)r, %(length)r)" % self.__dict__

class ReverseSequence(Sequence):
  """Complements of defined sequences"""
  def __init__(self, wc):
    self.name = wc.name + "*"
    self.prefix = wc.prefix
    self.full_name = wc.full_name + "*"
    self.length = wc.length
    self.seq = None # Stores the sequence once it has been defined.
    self.reversed = True
    self.dummy = wc.dummy
    self.wc = wc
  
  def __repr__(self):
    return "~" + repr(self.wc)
  
  def fix_seq(self, fixed_seq):
    """Constrian ourselves to a specific sequence."""
    self.wc.fix_seq(wc(fixed_seq))
  

class AnonymousSequence(Sequence):
  """Sequences we didn't lable and are thus anonymous."""
  num = 0
  def __init__(self, prefix, const, length=None):
    name = "_Anon"+repr(AnonymousSequence.num)
    Sequence.__init__(self, name, prefix, const, length)
    AnonymousSequence.num += 1


class SuperSequence(object):
  """Logical grouping of sequences"""
  def __init__(self, name, prefix, constraints, length=None):
    self.name = name
    self.prefix = prefix
    self.full_name = prefix + name
    self.seqs = []
    self.seq = None # Stores the sequence once it has been defined.
    self.reversed = False
    self.length = 0
    self.base_seqs = []
    wildcard = None
    # Process constraints
    for item in constraints:
      if isinstance(item, Sequence):
        # If item is a previously defined sequence, leave it alone, it's good
        self.seqs.append(item)
        self.base_seqs.append(item)
        self.length += item.length
      elif isinstance(item, SuperSequence):
        # Expand previous super-sequence
        self.seqs.append(item)
        self.base_seqs += item.base_seqs
        self.length += item.length
      else:
        # Otherwise it's a anonymous constraint
        try:
          anon_seq = AnonymousSequence(prefix, item)
          self.seqs.append(anon_seq)
          self.base_seqs.append(anon_seq)
          self.length += anon_seq.length
        except WildError:
          assert not wildcard, "Too many wildcards in super-sequence %s" % name
          wildcard = (len(self.seqs), len(self.base_seqs), item) # Index and entry of wildcard
    
    if not wildcard:
      assert self.length == length or length == None, "Length mismatch for sequence %s (%r != %r)" % (name, self.length, length)
    else:
      if length == None: raise WildError("Sequence %s has a ?. but no length specified" % name)
      wild_length = length - self.length  # Wildcard is set so that total length is right
      assert wild_length >= 0, "Sequence %s too short (%r > %r)" % (name, self.length, length)
      i, j, item = wildcard
      anon_seq = AnonymousSequence(prefix, item, wild_length)
      self.seqs.insert(i, anon_seq)
      self.base_seqs.insert(j, anon_seq)
      self.length += anon_seq.length
      assert self.length == length
    
    # Sequences of length 0 must be treated specially, they are not passed on 
    #   to lower levels, but are allowed for the convinience of users.
    if self.length == 0:
      self.dummy = True
    else:
      self.dummy = False
    
    self.wc = ReverseSuperSequence(self)
  
  def fix_seq(self, fixed_seq):
    """Constrian ourselves to a specific sequence."""
    assert len(fixed_seq) == self.length, \
            "Length of fixed sequence (%d) does not match length of %s (%d)" \
            % (len(fixed_seq), self.full_name, self.length)
    i = 0
    for seq in self.seqs:
      seq.fix_seq( fixed_seq[i:i+seq.length] )
      i += seq.length
  
  def __invert__(self):
    """Returns the Watson-Crick complementary sequence."""
    return self.wc
  def __repr__(self):
    return "SuperSequence(%(name)r, %(prefix)r, %(seqs)r, %(length)r)" % self.__dict__

class ReverseSuperSequence(SuperSequence):
  def __init__(self, wc):
    self.name = wc.name + "*"
    self.prefix = wc.prefix
    self.full_name = wc.full_name + "*"
    self.reversed = True
    self.length = wc.length
    self.seqs = [~seq for seq in wc.seqs[::-1]]
    self.seq = None # Stores the sequence once it has been defined.
    self.base_seqs = [~seq for seq in wc.base_seqs[::-1]]
    self.dummy = wc.dummy
    self.wc = wc

class Strand(SuperSequence):
  """Container for strands. Inherits from SuperSequence for convinience."""
  def __init__(self, name, prefix, constraints, length=None, dummy=False):
    SuperSequence.__init__(self, name, prefix, constraints, length)
    self.dummy = dummy
    
    # Keep track of whether this strand has been used in a structure yet.
    # If it never is, it may not be designed by some designers, so we warn user
    self.in_structure = False
    
  def __repr__(self):
    return "Strand(%(name)r, %(prefix)r, %(seqs)r, %(length)r, %(dummy)r)" % self.__dict__

class Structure(object):
  """Container for structures/complexes"""
  def __init__(self, name, prefix, strands, struct, opt=None):
    self.name = name
    self.prefix = prefix
    self.full_name = prefix + name
    self.opt = opt
    self.struct = struct
    self.mfe_struct = None # Stores the actual mfe structure once it's known
    self.strands = list(strands)
    self.seq = None # Stores the sequence once it has been defined.
    self.base_seqs = []
    sub_structs = [strand_struct for strand_struct in self.struct.split("+")] # Check that lengths match up
    assert len(strands) == len(sub_structs), "Mismatch: Structure %s is defined by %d strands, but secondary structure has %d strands" % (name, len(strands), len(sub_structs))
    for strand, sub_struct in zip(strands, sub_structs):
      assert isinstance(strand, Strand), "Structure %s must get strands" % name
      assert strand.length == len(sub_struct), "Mismatch: Strand %s in structure %s has length %d, but sub-structure %s implies %d" % (strand.name, name, strand.length, sub_struct, len(sub_struct))
      self.base_seqs += strand.base_seqs
  
  def fix_seq(self, fixed_seq):
    """Constrian ourselves to a specific sequence."""
    strand_seqs = fixed_seq.split("+")
    assert len(strand_seqs) == len(self.stands), "fix_seq: structure %s cannot have sequence fixed to %s because it has %d strands" % (self.name, fixed_seq, len(self.strands))
    for strand, strand_seq in zip(self.strands, strand_seqs):
      strand.fix_seq(strand_seq)
  
  def __repr__(self):
    return "Structure(%(name)r, %(prefix)r, %(strands)r, %(struct)r, %(opt)r)" % self.__dict__

class Kinetics(object):
  def __init__(self, name, prefix, inputs, outputs, low, high):
    self.name = name
    self.prefix = prefix
    self.full_name = prefix + name
    self.inputs = inputs
    self.outputs = outputs
    self.low = low
    self.high = high
  def __repr__(self):
    return "Kinetics(%(name)r, %(prefix)r, %(inputs)r, %(outputs)r, %(low)r, %(high)r)" % self.__dict__
