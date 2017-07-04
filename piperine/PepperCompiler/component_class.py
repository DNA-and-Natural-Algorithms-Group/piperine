import sys
import string

from utils import ordered_dict, PrintObject, error, warning
from DNA_classes import *
from component_parser import sequence_flag, domains_flag, nucleotide_flag

DEBUG = False

class Component(PrintObject):
  """Stores information for a DNA component"""
  def __init__(self, name, prefix, params):
    self.name = name
    self.prefix = prefix  # Prefix for all sequences, structs, etc. in this component.
    self.params = params
    
    self.seqs = ordered_dict()
    self.base_seqs = ordered_dict()  # Not super-sequences
    self.sup_seqs = ordered_dict()
    self.strands = ordered_dict()
    self.structs = ordered_dict()
    self.kinetics = ordered_dict()
    self.kin_num = 0
  
  def assert_(self, statement, message):
    """Raise error if statement is false. Adds component information to message."""
    prefix = "In component %s: " % self.name
    if not statement:
      error(prefix + str(message))
  
  ## Add information from document statements to object
  def add_sequence(self, name, const, length):
    if DEBUG: print "%s: sequence %s" % (self.name, name)
    self.assert_( name not in self.seqs, "Duplicate sequence definition for '%s'" % name )
    try:
      seq = Sequence(name, self.prefix, const, length)
    except AssertionError, e:
      self.assert_(False, str(e))
    self.base_seqs[name] = self.seqs[name] = seq
  
  def clean_const(self, old_const, name):
    """Replace all refferences to sequences with the actual sequences, expand domains, etc."""
    const = []
    for item in old_const:
      if item[0] == sequence_flag:
        seq_name, wc = item[1]
        self.assert_( seq_name in self.seqs, "Sequence '%s' referenced before definion (in sequence/strand '%s')" % (seq_name, name) )
        if not wc:
          seq = self.seqs[seq_name]
        else:
          seq = self.seqs[seq_name].wc
        const.append(seq)
      
      elif item[0] == domains_flag:
        seq_name, wc = item[1]
        self.assert_( seq_name in self.sup_seqs, "Sequence '%s' referenced before definion (in sequence/strand '%s')" % (seq_name, name) )
        if not wc:
          seq = self.sup_seqs[seq_name]
        else:
          seq = self.sup_seqs[seq_name].wc
        const += seq.seqs
      
      else:
        assert item[0] == nucleotide_flag, item
        const.append( item[1] )
    
    return const
  
  def add_super_sequence(self, name, const, length):
    if DEBUG: print "%s: super-sequence %s" % (self.name, name)
    self.assert_( name not in self.seqs, "Duplicate sequence definition for '%s'" % name )
    const = self.clean_const(const, name)
    try:
      seq = SuperSequence(name, self.prefix, const, length)
    except AssertionError, e:
      self.assert_(False, str(e))
    self.sup_seqs[name] = self.seqs[name] = seq
    # Add anonymous sequences
    for seq in self.sup_seqs[name].seqs:
      if seq.reversed: # Look for the standard form of the sequences
        seq = seq.wc
      if seq.name not in self.seqs:
        self.assert_(isinstance(seq, AnonymousSequence), "In super-sequence %s, sequence %s has not been defined yet." % (name, seq.name))
        self.seqs[seq.name] = seq
        self.base_seqs[seq.name] = seq
  
  def add_strand(self, dummy, name, const, length):
    if DEBUG: print "%s: strand %s" % (self.name, name)
    self.assert_( name not in self.strands, "Duplicate strand definition for '%s'" % name )
    const = self.clean_const(const, name)
    try:
      self.strands[name] = strand = Strand(name, self.prefix, const, length, dummy)
    except AssertionError, e:
      self.assert_(False, str(e))
    self.assert_(strand.length > 0, "Strand %s was defined with length 0" % name)
    # Add anonymous sequences we just created to the lists
    for seq in strand.seqs:
      if seq.reversed: # Look for the standard form of the sequences
        seq = seq.wc
      if seq.name not in self.seqs:
        self.assert_(isinstance(seq, AnonymousSequence), "In strand %s, sequence %s has not been defined yet." % (name, seq.name))
        self.seqs[seq.name] = seq
        self.base_seqs[seq.name] = seq
    for seq in strand.base_seqs:
      if seq.reversed:
        seq.wc.in_strand = True
      else:
        seq.in_strand = True
  
  def add_structure(self, opt, name, strands, struct):
    if DEBUG: print "%s: structure %s" % (self.name, name)
    self.assert_( name not in self.structs, "Duplicate structure definition for '%s'" % name )
    
    # Convert from list of strand names to list of strands
    for n, strand in enumerate(strands):
      self.assert_( strand in self.strands, "Strand '%s' referenced before definion (in structure '%s')" % (strand, name) )
      strands[n] = self.strands[strand]
      strands[n].in_structure = True
    
    isdomain, struct = struct
    if isdomain: # This is a domain-based structure
      sub_structs = struct.split("+")
      full_struct = ""
      self.assert_( len(sub_structs) == len(strands), "Mismatched number of strands: Structure %s has %d strands, but structures %s implies %d." % (name, len(strands), struct, len(sub_structs)) )
      # For each strand expand out the structure
      for sub_struct, strand in zip(sub_structs, strands):
        self.assert_( len(sub_struct) == len(strand.seqs), "Mismatch: Strand %s in structure %s has %d domain(s), but sub-structure %s implies %d" % (strand.name, name, len(strand.seqs), sub_struct, len(sub_struct)) )
        for dp, domain in zip(sub_struct, strand.seqs):
          full_struct += dp * domain.length
        full_struct += "+"
      struct = full_struct[:-1] # Get rid of trailing +
    try:
      self.structs[name] = Structure(name, self.prefix, strands, struct, opt)
    except AssertionError, e:
      self.assert_(False, str(e))
  
  def add_kinetic(self, low, high, inputs, outputs):
    if DEBUG: print "%s: kinetic Kin%d" % (self.name, self.kin_num)
    for n, struct in enumerate(inputs):
      self.assert_( struct in self.structs, "Kinetic statement uses structure '%s' before it is defined." % struct )
      inputs[n] = self.structs[struct]
    for n, struct in enumerate(outputs):
      self.assert_( struct in self.structs, "Kinetic statement uses structure '%s' before it is defined." % struct )
      outputs[n] = self.structs[struct]
    if not low:
      low = 0
    if not high:
      high = float("inf")
    
    name = "Kin%d" % self.kin_num
    self.kin_num += 1
    self.kinetics[name] = Kinetics(name, self.prefix, list(inputs), list(outputs), low, high)
  
  def add_IO(self, inputs, outputs):
    """Add I/O information once we've read the component."""
    self.input_seqs = []
    self.input_structs = []
    for (seq_name, wc), struct_name in inputs:
      self.assert_( seq_name in self.seqs, "Declare statement references undefined sequence '%s'" % seq_name )
      if wc:
        self.input_seqs.append( self.seqs[seq_name].wc )
      else:
        self.input_seqs.append( self.seqs[seq_name] )
      
      if struct_name:
        self.assert_( struct_name in self.structs, "Declare statement references undefined structure '%s'" %  struct_name )
        self.input_structs.append( self.structs[struct_name] )
      else:
        self.input_structs.append(None)
    
    self.output_seqs = []
    self.output_structs = []
    for (seq_name, wc), struct_name in outputs:
      self.assert_( seq_name in self.seqs, "Declare statement references undefined sequence '%s'" % seq_name )
      if wc:
        self.output_seqs.append( self.seqs[seq_name].wc )
      else:
        self.output_seqs.append( self.seqs[seq_name] )
      
      if struct_name:
        self.assert_( struct_name in self.structs, "Declare statement references undefined structure '%s'" %  struct_name )
        self.output_structs.append( self.structs[struct_name] )
      else:
        self.output_structs.append(None)
        
  
  
  ## Outputs
  def output_synthesis(self, prefix, outfile):
    """Output synthesis of all data into a single file."""
    if prefix:
      outfile.write("#\n## Component %s\n" % prefix[:-1])
    else:
      outfile.write("#\n## Top Component\n")
    
    # Define sequences
    for seq in self.base_seqs.values():
      if not seq.dummy:
        if not seq.in_strand:
          warning("Sequence %s is defined but never used in a strand. It probably will not be designed." % seq.full_name)
        outfile.write("sequence %s = %s : %d\n" % (seq.full_name, seq.const, seq.length))
    
    # Define super-sequences
    for sup_seq in self.sup_seqs.values():
      if not seq.dummy:
        const = string.join([seq.full_name for seq in sup_seq.seqs if not seq.dummy], " ")
        outfile.write("sup-sequence %s = %s : %d\n" % (sup_seq.full_name, const, sup_seq.length))
    
    # Define strands
    for strand in self.strands.values():
      if not strand.in_structure:
        warning("Strand %s is defined but never used in a structure. It may not be designed." % strand.full_name)
      const = string.join([seq.full_name for seq in strand.seqs if not seq.dummy], " ")
      if strand.dummy:
        dummy = "[dummy] "
      else:
        dummy = ""
      outfile.write("strand %s%s = %s : %d\n" % (dummy, strand.full_name, const, strand.length))
    
    # Define structures
    for struct in self.structs.values():
      strands = string.join([strand.full_name for strand in struct.strands], " + ")
      outfile.write("structure [%dnt] %s = %s : %s\n" % (struct.opt, struct.full_name, strands, struct.struct))
    
    # Define kinetics
    for kin in self.kinetics.values():
      inputs = string.join([struct.full_name for struct in kin.inputs], " + ")
      outputs = string.join([struct.full_name for struct in kin.outputs], " + ")
      outfile.write("kinetic [%f /M/s < k < %f /M/s] %s -> %s\n" % (kin.low, kin.high, inputs, outputs))
  
  
  def output_nupack(self, prefix, outfile):
    """Compile data into NUPACK format and output it"""
    if prefix:
      outfile.write("#\n## Component %s\n" % prefix[:-1])
    else:
      outfile.write("#\n## Top Component\n")
    
    used_seqs = set()
    
    # Define structures
    for struct in self.structs.values():
      outfile.write("structure %s = %s\n" % (struct.full_name, struct.struct))
      # Add all sequences in this structure to set of used sequences.
      used_seqs.update([ x for x in struct.base_seqs if not x.reversed] + \
                       [~x for x in struct.base_seqs if x.reversed])
  
    # Define sequences
    for seq in self.base_seqs.values():
      self.assert_(isinstance(seq, Sequence), "Expected Sequence object instead of %r" % seq)
      if seq not in used_seqs:
        warning("Sequence %s is defined, but never used in a structure. It may not be designed." % seq.full_name)
      
      if not seq.dummy:
        outfile.write("sequence %s = %s\n" % (seq.full_name, seq.const))
    
    # Apply sequences to structures and set objective function
    for struct in self.structs.values():
      seqs = string.join([seq.full_name for seq in struct.base_seqs if not seq.dummy])
      outfile.write("%s : %s\n" % (struct.full_name, seqs))
      if struct.opt: # Optimization parameter
        outfile.write("%s < %f\n" % (struct.full_name, struct.opt))
