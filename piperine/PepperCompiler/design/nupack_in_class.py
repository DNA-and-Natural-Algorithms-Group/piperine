from DNA_nupack_classes import *

from ..utils import ordered_dict, PrintObject

class Spec(PrintObject):
  def __init__(self):
    self.structs = ordered_dict()
    self.seqs = ordered_dict()
  
  def add_structure(self, name, struct):
    assert name not in self.structs, "Duplicate structure definition"
    self.structs[name] = Structure(name, struct)
  def add_sequence(self, name, const):
    assert name not in self.seqs, "Duplicate sequence definition"
    num = len(self.seqs) # Number the sequences consecutively.
    self.seqs[name] = Sequence(name, const, num)
  def add_apply(self, struct_name, seq_names):
    #print self.structs
    struct = self.structs[struct_name]
    seqs = [self.load_seq(seq_name) for seq_name in seq_names]
    struct.set_seqs(seqs)
  
  def load_seq(self, foo):
    name, parity = foo
    if parity == "":
      return self.seqs[name]
    else: # parity == "*"
      return ~self.seqs[name]
