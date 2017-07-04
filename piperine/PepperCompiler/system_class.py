"""
The System class stores all of the information in a system file.
"""

import string

import DNA_classes
from utils import ordered_dict, default_ordered_dict, PrintObject, error

DEBUG = False

def load_file(basename, args, prefix, path=".", includes=None):
  """Load the file basename.(sys|comp) with args from path."""
  # Imported in the function to avoid circular import error.
  import os
  from system_parser import load_system
  from component_parser import load_component
 
  if includes:
    paths = [path] + includes
  else:
    paths = [path]
  
  for incpath in paths:
    basenameandpath = os.path.join(incpath, basename) # Add the correct directory name
    new_path = os.path.dirname(basenameandpath) # Local directory of basename
    
    sys_name = basenameandpath+".sys"   # Name if file is a system specification
    comp_name = basenameandpath+".comp" # Name if file is a component spec.
    
    issys = os.path.isfile(sys_name)   # Is it a system file?
    iscomp = os.path.isfile(comp_name) # Is it a component file?

    if issys and iscomp:
      error("Ambiguous specification: Both '%s' and '%s' exist. Please remove the one that does not belong and rerun the compiler." % (sys_name, comp_name))

    if issys or iscomp:
      break
  
  # Check that exactly one of the two types exists
  if not (issys or iscomp):
    error("Neither '%s' nor '%s' exist" % (sys_name, comp_name))
  
  # And load it
  if issys:
    return load_system(sys_name, args, prefix, new_path, includes=includes)
  
  else: # iscomp
    return load_component(comp_name, args, prefix)

class System(PrintObject):
  """Stores all the information in a system's connectivity file"""
  def __init__(self, path, name, prefix, params, includes=None):
    """Initialized the cicuit with the declare statement"""
    self.path = path  # Filesystem path to load components relative to.
    self.name = name
    self.prefix = prefix  # Prefix for all sub-components.
    self.includes = includes
    
    self.template = ordered_dict()
    self.signals = ordered_dict()
    self.lengths = ordered_dict()
    self.signal_structs = default_ordered_dict(list, call=True)  # Structures output as signals
    self.signal_dummy_structs = default_ordered_dict(list, call=True)  # Dummy structures marked as inputs of signals
    self.components = ordered_dict()
    
    # Pointers to subcomponent's objects
    self.seqs = ordered_dict()
    self.base_seqs = ordered_dict()  # Not super-sequences
    self.sup_seqs = ordered_dict()
    self.strands = ordered_dict()
    self.structs = ordered_dict()
    self.kinetics = ordered_dict()
  
  ## Add information from document statements to object
  def add_import(self, imports):
    if DEBUG: print "import", imports
    for path, name in imports:
      if name == None:
        # filename is used as the internal name by default
        if "/" not in path:
          name = path
        else:
          name = path.split("/")[-1]  # Strip off lower directories
      if DEBUG: print "import", name, path
      assert name not in self.template, "Duplicate import %s" % name
      self.template[name] = path

  def add_component(self, comp_name, templ_name, templ_args, inputs, outputs):
    if DEBUG: print "component", comp_name, templ_name, templ_args, inputs, outputs
    # Setup components
    assert templ_name in self.template, "Template referenced before import: " + templ_name
    assert comp_name not in self.components, "Duplicate component definition: " + comp_name
    self.components[comp_name] = this_comp = load_file(self.template[templ_name], templ_args, self.prefix + comp_name + "-", self.path, includes=self.includes)
    assert len(inputs) == len(this_comp.input_seqs),   "Length mismatch. %s / %s: %r != %r" % (comp_name, templ_name, len( inputs), len(this_comp.input_seqs ))
    assert len(outputs) == len(this_comp.output_seqs), "Length mismatch. %s / %s: %r != %r" % (comp_name, templ_name, len(outputs), len(this_comp.output_seqs))
    # Constrain all component inputs and outputs
    ### TODO: marry these 2 together in a more eligent way.
    if isinstance(this_comp, System): # If it's actually a system
      for (glob_name, glob_wc), (loc_name, loc_wc) in zip(list(inputs)+list(outputs), this_comp.input_seqs+this_comp.output_seqs):
        wc = (glob_wc != loc_wc)  # Are these signals complementary?
        if glob_name not in self.signals:
          self.signals[glob_name] = [(loc_name, comp_name, wc)]
          self.lengths[glob_name] = this_comp.lengths[loc_name]
        else:
          self.signals[glob_name].append( (loc_name, comp_name, wc) )
          assert self.lengths[glob_name] == this_comp.lengths[loc_name]
      # Collect structures that could represent signals
      for (glob_name, glob_wc), loc_structs in zip(list(outputs), this_comp.output_structs):
        # TODO: deal with the wc aspect of this appropriately.
        self.signal_structs[glob_name] += loc_structs
      # ... and point all dummy inputs to those actual structures.
      for (glob_name, glob_wc), loc_structs in zip(list(inputs), this_comp.input_structs):
        # TODO: deal with the wc aspect of this appropriately.
        self.signal_dummy_structs[glob_name] += loc_structs
        for loc_struct in loc_structs:
          loc_struct.actual_structs = self.signal_structs[glob_name]
    else: # Otherwise it's a component, so we want to constrain sequences
      for (glob_name, glob_wc), loc_seq in zip(list(inputs)+list(outputs), this_comp.input_seqs+this_comp.output_seqs):
        wc = (glob_wc != loc_seq.reversed)  # Are these signals complementary?
        if glob_name not in self.signals:
          self.signals[glob_name] = [(loc_seq, comp_name, wc)]
          assert not loc_seq.dummy, "In system %s: Signal %s represented by a dummy (length 0) sequence %s. Dummy signals are not allowed." % (self.name, glob_name, loc_seq.name)
          self.lengths[glob_name] = loc_seq.length
        else:
          self.signals[glob_name].append( (loc_seq, comp_name, wc) )
          assert self.lengths[glob_name] == loc_seq.length, "In system %s: Lengths of signal %s is inconsistent between %s and %s (%d != %d)" % (self.name, glob_name, self.signals[glob_name][0][0].full_name, loc_seq.full_name, self.lengths[glob_name], loc_seq.length)
      # Collect structures that could represent signals
      for (glob_name, glob_wc), loc_struct in zip(list(outputs), this_comp.output_structs):
        # TODO: deal with the wc aspect of this appropriately.
        self.signal_structs[glob_name].append(loc_struct)
      # ... and point all dummy inputs to those actual structures.
      for (glob_name, glob_wc), loc_struct in zip(list(inputs), this_comp.input_structs):
        # TODO: deal with the wc aspect of this appropriately.
        if loc_struct:
          self.signal_dummy_structs[glob_name].append(loc_struct)
          loc_struct.actual_structs = self.signal_structs[glob_name]
    
    # Point to all objects in the component
    # For each type of object: seqs, strands, ...
    for type_ in "seqs", "base_seqs", "sup_seqs", "strands", "structs", "kinetics":
      comp_objs = this_comp.__dict__[type_] # this_comp.seqs, this_comp.base_seqs, ...
      system_objs = self.__dict__[type_]   # self.seqs, self.base_seqs, ...
      # Point to all of those items from here with a prefix added to the name
      for name, obj in comp_objs.items():
        system_objs[comp_name + "-" + name] = obj
  
  def add_IO(self, inputs, outputs):
    """Add I/O information once we've read the component."""
    self.input_seqs = []
    self.input_structs = []
    for seq_name, wc in inputs:
      assert seq_name in self.signals
      self.input_seqs.append( (seq_name, wc) )
      self.input_structs.append( self.signal_dummy_structs[seq_name] )
    
    self.output_seqs = []
    self.output_structs = []
    for seq_name, wc in outputs:
      assert seq_name in self.signals, seq_name + str(self.signals)
      self.output_seqs.append( (seq_name, wc) )
      self.output_structs.append( self.signal_structs[seq_name] )
  
  
  def output_synthesis(self, prefix, outfile):
    """Output synthesis of all data into a single file."""
    if prefix:
      outfile.write("#\n## Subsystem %s\n" % prefix[:-1])
    else:
      outfile.write("#\n## Top System\n")
    # For each component write it's contents with prefix "name-".
    for component_name, template in self.components.items():
      template.output_synthesis(prefix+component_name+"-", outfile)
    
    # For each signal sequence connecting components constrain them to be be equal.
    if prefix:
      outfile.write("#\n## System %s Connectors\n" % prefix[:-1])
    else: 
      outfile.write("#\n## Top System Connectors\n")
    for signal in self.signals:
      length = self.lengths[signal]
      signal_name = prefix + signal
      # Create signal sequence object
      outfile.write("#\n## Signal %s\n" % signal_name)
      outfile.write("sequence %s = %s : %d\n" % (signal_name, "N" * length, length))
      
      # Constrain all sequences representing that signal
      outfile.write("equal %s " % signal_name)
      for loc_seq, component_name, wc in self.signals[signal]:
        if isinstance(loc_seq, (DNA_classes.Sequence, DNA_classes.SuperSequence)):
          loc_name = loc_seq.full_name
        else: # it's a system signal, so it's just a name
          loc_name = prefix + component_name + "-" + loc_seq
        if wc:
          outfile.write("%s* " % loc_name)
        else:
          outfile.write("%s " % loc_name)
      outfile.write("\n")
  
  def output_nupack(self, prefix, outfile):
    """Compile data into NUPACK format and output it"""
    if prefix:
      outfile.write("#\n## Subsystem %s\n" % prefix[:-1])
    else:
      outfile.write("#\n## Top System\n")
    # For each component write it's contents in Zadeh's format with prefix "name-".
    for comp_name, template in self.components.items():
      template.output_nupack(prefix+comp_name+"-", outfile)
    
    # For each signal sequence connecting components constrain them to be be equal.
    # To force this constraint I make them all complementary to a single dummy strand
    if prefix:
      outfile.write("#\n## System %s Connectors\n" % prefix[:-1])
    else: 
      outfile.write("#\n## Top System Connectors\n")
    for signal in self.signals:
      length = self.lengths[signal]
      signal_name = prefix + signal
      wc_name = signal_name + "-_WC"  # A dummy variable wc complement to signal
      
      outfile.write("#\n## Signal %s\n" % signal_name)
      outfile.write("sequence %s = %s\n" % (signal_name, "N" * length))
      outfile.write("sequence %s = %s\n" % (wc_name, "N" * length))
      
      dummy_name = "%s-_Self" % signal_name
      outfile.write("structure %s = %s\n" % (dummy_name, "(" * length + "+" + ")" * length))
      outfile.write("%s : %s %s\n" % (dummy_name, wc_name, signal_name))
      
      # For each instance of the signal sequence, build a structure to 
      #   constrain it to the master signal sequence
      for loc_seq, comp_name, wc in self.signals[signal]:
        if isinstance(loc_seq, DNA_classes.Sequence):
          sig_name = comp_name + "-" + loc_seq.name
          seqs = prefix + sig_name
        elif isinstance(loc_seq, DNA_classes.SuperSequence):
          # If it's a super-seq, list the subsequences
          sig_name = comp_name + "-" + loc_seq.name
          seqs = string.join([prefix+comp_name+"-"+seq.name for seq in loc_seq.base_seqs if not seq.dummy])
        else: # it's a system signal
          sig_name = comp_name + "-" + loc_seq
          seqs = prefix + sig_name
        
        dummy_name = signal_name + "-" + sig_name
        outfile.write("structure %s = %s\n" % (dummy_name, "(" * length + "+" + ")" * length))
        
        if wc: # If it's complementary to the signal, then we can enforce that directly
          outfile.write("%s : %s %s\n" % (dummy_name, signal_name, seqs))
        else:  # If it's equal, then we must force it complementary to the complement 'wc_name'
          outfile.write("%s : %s %s\n" % (dummy_name, wc_name, seqs))
