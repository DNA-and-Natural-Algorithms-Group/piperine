from PIL_DNA_classes import *

from ..utils import ordered_dict, PrintObject

def get_seqs(seq_names, seq_dict):
  """Get list of sequences from list of names."""
  seqs = []
  for name in seq_names:
    if name[-1] == "*": # If it's a complement of a sequence
      name = name[:-1]
      assert name in seq_dict, "Sequence %s referenced before definition" % name # Make sure it's really a sequence
      seqs.append(seq_dict[name].wc)
    else:
      assert name in seq_dict, "Sequence %s referenced before definition" % name
      seqs.append(seq_dict[name])
  return seqs

def get_strands(names, strand_dict):
  """Get list of sequences from list of names."""
  strands = []
  for name in names:
    assert name in strand_dict, "Strand %s referenced before definition" % name
    strands.append(strand_dict[name])
  return strands

class Spec(PrintObject):
  def __init__(self):
    self.base_seqs = ordered_dict()
    self.sup_seqs = ordered_dict()
    self.seqs = ordered_dict() # All sequences
    self.strands = ordered_dict()
    self.structs = ordered_dict()
    self.equals = []
  
  def add_seq(self, name, template):
    assert name not in self.seqs, "Duplicate sequence definition: %s" % name
    self.base_seqs[name] = Sequence(name, template)
    self.seqs[name] = self.base_seqs[name]
  
  def add_sup_seq(self, name, sub_seq_names):
    assert name not in self.seqs, "Duplicate sequence definition: %s" % name
    sub_seqs = get_seqs(sub_seq_names, self.seqs)
    self.sup_seqs[name] = SuperSequence(name, sub_seqs)
    self.seqs[name] = self.sup_seqs[name]
  
  def add_strand(self, name, seq_names, dummy):
    assert name not in self.strands, "Duplicate strand definition: %s" % name
    seqs = get_seqs(seq_names, self.seqs)
    self.strands[name] = Strand(name, seqs, dummy)
  
  def add_struct(self, name, strand_names, struct, params):
    assert name not in self.structs, "Duplicate structure definition: %s - %r" % (name, (strand_names, struct, params))
    strands = get_strands(strand_names, self.strands)
    self.structs[name] = Structure(name, strands, struct, params)
  
  def add_equal(self, seq_names):
    seqs = get_seqs(seq_names, self.seqs)
    l = seqs[0].length
    for seq in seqs:
      assert seq.length == l, "equal statement has different length sequences: %r & %r" % (seqs[0].name, seq.name)
    self.equals.append( seqs )

