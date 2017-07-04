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
  def __init__(self, name, template_seq, num):
    self.name = name
    self.seq = template_seq
    self.nseq = None
    self.num = num
    self.reversed = False
    self.length = len(self.seq)
    # Build the dummy sequence for the W-C complement
    self.wc = ReverseSequence(self)
  def __invert__(self):
    """Returns the Watson-Crick complementary sequence."""
    return self.wc

class ReverseSequence(Sequence):
  """Complements of defined sequences"""
  def __init__(self, wc):
    self.name = wc.name + "*"
    self.seq = seq_comp(wc.seq)
    self.nseq = None
    self.length = wc.length
    self.num = wc.num
    self.reversed = True
    self.wc = wc

def get_bonds(struct):
  bonds = []
  stack = []
  i = 0
  for symb in struct:
    if symb == "(":
      stack.append(i)
    if symb == ")":
      assert len(stack) != 0
      start = stack.pop() # The most recent open-paren
      bonds.append( (start, i) )
    
    # 'i' is the index (but index doesn't include strand breaks)
    if symb != "+":
      i += 1
  
  return bonds

class Structure(object):
  """Container for structures/complexes"""
  def __init__(self, name, struct):
    self.name = name
    self.struct = struct
    self.length = len(struct.replace("+", ""))
    self.bonds = get_bonds(struct)
    self.seqs = None
    self.seq = None
    self.nseq = None
  def set_seqs(self, seqs):
    """Set the sequences for a structure."""
    assert not self.seqs, "Sequences have already been set for this structure."
    self.seqs = tuple(seqs)
    assert self.length == sum([seq.length for seq in self.seqs])
    self.seq = string.join([seq.seq for seq in self.seqs], "")

  def seq_loc(self, index):
    """Find out which sequence of the structure that the index falls into."""
    num = 0  # Current seq number
    #print self.name, index
    while self.seqs[num].length  <=  index:
      index -= self.seqs[num].length  # Move past sequence
      num += 1

    if not self.seqs[num].reversed: # If it's normal direction
      return self.seqs[num].num, index, True
    else:  # It's backwards
      return self.seqs[num].num, self.seqs[num].length - index - 1, False

  def __repr__(self):
    return "Structure(%(name)r, %(struct)r)" % self.__dict__
