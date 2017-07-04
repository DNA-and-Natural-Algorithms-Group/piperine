from __future__ import division

import os

from utils import ordered_dict, PrintObject, error
import nupack_out_grammar as ngram
from multistrand import DNAkinfold

def read_design(filename):
  """Extracts the designed sequences and the mfe structures"""
  if not os.path.isfile(filename):
    error("Cannot load design. No such file '%s'." % filename)
  stats, total_n_star = ngram.document.parseFile(filename)
  seqs = {}
  structs = {}
  # Catches everything (struct, seq and seq*) except bad inters (struct N struct)
  #print stats
  for stat in stats:
    #print stat
    name, seq, n_star, gc, mfe, ideal_struct, actual_struct = stat
    seqs[name] = seq
    structs[name] = actual_struct
  return seqs

class Complex(PrintObject):
  def __init__(self, strands, struct):
    self.strands = strands; self.struct = struct

def test_kinetics(kin, cleanup, trials=24, time=100000, temp=25, conc=1.0, out_interval=-1):
  """Test times for inputs to combine/separate into outputs"""
  used_strands = ordered_dict()
  ins = []
  for struct in kin.inputs:
    # And add the starting complexes/structures
    strand_names = [strand.full_name for strand in struct.strands]
    ins.append( Complex(strand_names, struct.struct) )
    # Keep track of strands that will be used
    for strand in struct.strands:
      assert strand.full_name not in used_strands, (strand, used_strands) # Sanity check
      used_strands[strand.full_name] = strand.seq
  
  outs = []
  for struct in kin.outputs:
    strand_names = [strand.full_name for strand in struct.strands]
    # TODO: check for conservation of strands
    outs.append( Complex(strand_names, "DISASSOC") )
  # HACK: Multistrand doesn't currently work with multiple DISASSOC structures, 
  #   we just wait for the first one to form.
  outs = [ outs[0] ]
  # HACK: Likewise the backwords reaction must only have one DISASSOC.
  back = [ Complex(ins[0].strands, "DISASSOC") ]
  
  return DNAkinfold(used_strands, ins, back, outs, trials, time, temp, conc, out_interval, cleanup)

def test_spuradic(structs, cleanup, trials=24, time=10, temp=25, conc=1.0, out_interval=-1):
  """Test kinetics for spuratic interaction between structures."""
  used_strands = ordered_dict()
  ins = []
  for i, struct in enumerate(structs):
    # We must be able to distinguish two strands which are actually the same 
    #   strand and thus have the same name. Thus we add a struct # i at the end.
    i = str(i)
    # And add the starting complexes/structures
    strand_names = [strand.full_name + i for strand in struct.strands]
    ins.append( Complex(strand_names, struct.struct) )
    # Keep track of strands that will be used
    for strand in struct.strands:
      name = strand.full_name + i
      assert name not in used_strands
      used_strands[name] = strand.seq
  
  outs = [] # We are not testing for any forward reactions.
  # HACK: Likewise the backwords reaction must only have one DISASSOC.
  back = [Complex(ins[0].strands, "DISASSOC")]
  
  return DNAkinfold(used_strands, ins, back, outs, trials, time, temp, conc, out_interval, cleanup)
