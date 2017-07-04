#!/usr/bin/env python
from __future__ import division

import sys
from utils import error
if sys.version_info < (2, 5):
  error("Must use python 2.5 or greater.")

import string
from copy import copy
from subprocess import CalledProcessError

from compiler import load
from kinetics import read_design, test_kinetics, test_spuradic

from system_class import System
from component_class import Component
from DNA_classes import wc
import myStat as stat

def get_components(obj, prefix=""):
  """Return all components and the prefixes showing how to reach them."""
  if isinstance(obj, Component):
    return [(obj, prefix)]
  elif isinstance(obj, System):
    components = []
    for component_name, component in obj.components.items():
      components += get_components(component, prefix + component_name + "-")
    return components
  else:
    raise Exception, 'Object "%r" is neither a System or Component.' % obj

def finish(savename, designname, seqsname, strandsname, run_kin, cleanup, trials, time, temp, conc, spurious, spurious_time):
  """
  Finish compiling a specification.
  
  Currently, it
    1) produces a .seqs file of all designed sequences,
    2) optionally produces a "strands to order" file and
    3) tests the specified kinetic paths
  
  In the future it might test for bad kinetics, etc.
  """
  
  print "Finishing compilation of %s ..." % savename
  # Re-load the design system/component
  system = load(savename)
  
  print "Applying the design from '%s'" % designname
  seqs = read_design(designname)
  apply_design(system, seqs)
  
  # Document all sequences, super-sequences, strands, and structures
  print "Writing sequences file: %s" % seqsname
  f = open(seqsname, "w")
  
  f.write("# Sequences\n")
  for name, seq in system.seqs.items():
    f.write("sequence %s = %s\n" % (name, seq.seq))
  f.write("# Strands\n")
  for name, strand in system.strands.items():
    f.write("strand %s = %s\n" % (name, strand.seq))
  f.write("# Structures\n")
  for name, struct in system.structs.items():
    f.write("structure %s = %s\n" % (name, struct.seq))
  f.close()
  
  # Document all strands that will be used in the final experiment
  if strandsname:
    print 'Writing "stands to order" file: %s' % strandsname
    f = open(strandsname, "w")
    for name, strand in system.strands.items():
      if not strand.dummy:
        f.write("strand %s\t%s\n" % (name, strand.seq) )
    f.close()
  
  # TODO: Write a thermodynamic scorecard.
  
  # Run Kinetic tests
  if run_kin:
    print "Testing Kinetics."
    print "Trials: %d" % trials
    print "SimTime: %.1f s" % time
    print "Temperature: %.1f deg C" % temp
    print "Concentration: %f uM" % conc
    for component, prefix in get_components(system):
      kinetic(component, prefix, cleanup, trials, time, temp, conc)
  
  # Test spuradic kinetics
  if spurious:
    print "Input structures to be tested for spuradic kinetics."
    for struct1 in system.structs.values():
      for struct2 in system.structs.values():
        print "Testing spurious kinetics of:", struct1.full_name, struct2.full_name
        try:
          process_kinetics(test_spuradic([struct1, struct2], cleanup, trials, spurious_time, temp, conc))
        except CalledProcessError, e:
          error(str(e))


def apply_design(system, seqs):
  """Assigns designed sequences and provided mfe structures to the respective objects."""
  # Assign all the designed sequences
  for name, seq in system.base_seqs.items():
    if not seq.dummy:
      seq.seq  = seqs[name]
      assert len(seq.seq) == seq.length
      seq.wc.seq = wc(seq.seq)
      assert seq.wc.seq == seqs[name + "*"]
    else:
      seq.seq = seq.wc.seq = ""
  
  for sup_seq in system.sup_seqs.values():
    sup_seq.seq  = string.join([seq.seq for seq in sup_seq.base_seqs], "")
    sup_seq.wc.seq = string.join([seq.seq for seq in sup_seq.wc.base_seqs], "")
  
  for strand in system.strands.values():
    strand.seq = string.join([seq.seq for seq in strand.base_seqs], "")
  
  for name, struct in system.structs.items():
    struct.seq = string.join([strand.seq for strand in struct.strands], "+")
    assert struct.seq == seqs[name], "Design is inconsistant! %s != %s" % (struct.seq, seqs[name])


def choices(x):
  """
  Takes a list of lists x and returns all possible choices of one element from each sublist.
  
  >>> choices([range(2), range(3), ["foo"]])
  [ [0, 0, "foo"], [0, 1, "foo"], [0, 2, "foo"], [1, 0, "foo"], [1, 1, "foo"], [1, 2, "foo"] ]
  """
  if len(x) == 0:
    return [[]]
  
  else:
    res = []
    for item in x[0]:
      res += [ [item] + sub for sub in choices(x[1:]) ]
    return res

# TODO: get rid of need for component, so get rid of recursion.
def kinetic(component, prefix, cleanup, trials, time, temp, conc):
  """Test all kinetic pathways in a component."""
  for kin in component.kinetics.values():
    ## First, we gather the set of all possible input structures.
    pos_inputs = [None] * len(kin.inputs)
    for n, struct in enumerate(kin.inputs):
      # If this is a dummy input structure and there are actual structures that will replace it
      #   AND: if the dummy structure is just a single strand (otherwise, we may not know what to do).
      if struct in component.input_structs and struct.actual_structs and len(struct.strands) == 1:
        # ... then we list these structures
        pos_inputs[n] = struct.actual_structs
      # Otherwise, we just use the dummy default structure.
      else:
        pos_inputs[n] = [struct]
    
    ## Next, we test each permutation of input structures.
    print_kin(kin, prefix[:-1])
    for inputs in choices(pos_inputs):
      new_kin = sub_inputs(kin, inputs)
      print_kin(new_kin, prefix[:-1])
      res = test_kinetics(new_kin, cleanup, trials, time, temp, conc)
      process_kinetics(res)

def sub_inputs(kin, inputs):
  """Substitute actual strands from structures in 'inputs' for the dummy strands in kin"""
  new_kin = copy(kin)
  ## Replacing inputs is easy
  new_kin.inputs = inputs
  
  ## Replacing outputs requires searching through and substituting strands
  new_kin.outputs = copy(kin.outputs)
  # First create a dictionary of all replacements
  replace = {}
  assert len(inputs) == len(kin.inputs)
  for real_struct, dummy_struct in zip(inputs, kin.inputs):
    if dummy_struct != real_struct:
      assert len(dummy_struct.strands) == 1, repr(dummy_struct)
      dummy_strand = dummy_struct.strands[0]
      replace[dummy_strand.full_name] = real_struct.strands
  
  # Then go through the outputs and actually replace them
  for n, struct in enumerate(new_kin.outputs):
    new_struct = copy(struct)
    new_kin.outputs[n] = new_struct
    new_struct.struct = None  # HACK: We don't use the structure for now, so we don't have to deal with this.
    new_struct.strands = []
    for strand in struct.strands:
      if strand.full_name in replace:
        new_struct.strands += replace[strand.full_name]
      else:
        new_struct.strands.append(strand)
  
  return new_kin

def print_kin(kin, gate_name):
  """Print kinetic testing info"""
  print
  print "kinetic", gate_name, ":",
  for struct in kin.inputs:
    print struct.full_name,
  print "->",
  for struct in kin.outputs:
    print struct.full_name,
  print
  sys.stdout.flush()

def process_kinetics(ret):
  """Process results of multistrand and print summary, etc."""
  forward, reverse, overtime, summary = ret
  # Process results
  num_for = len(forward)
  num_rev = len(reverse)
  num_over = len(overtime)
  num_trials = num_for + num_rev + num_over
  
  if num_over > 0:
    print
    print "WARNING: %d/%d trajectories went overtime." % (num_over, num_trials)
    print "Reported statistics may be unreliable."
    print
  
  # Rate at which collisions that will eventually go forward happen
  for_coll_rate = sum([coll_rate for (time, coll_rate) in forward]) / num_trials
  
  for_rates = [1/time for (time, coll_rate) in forward]
  for_rate = stat.mean(for_rates)
  for_rate_stddev = stat.stddev(for_rates)
  
  print "* %d/%d trajectories went forward." % (num_for, num_trials)
  if num_for > 0:
    print "  Estimated Forward Collision Rate: %f /uM/s" % (for_coll_rate / 1000000)
    print "  Estimated Forward Trajectory Rate: %f (std-dev %f) /s" % (for_rate, for_rate_stddev)
    print
  
  # Rate at which collisions that will eventually reverse happen
  rev_coll_rate = sum([coll_rate for (time, coll_rate) in reverse]) / num_trials
  
  rev_rates = [1/time for (time, coll_rate) in reverse]
  rev_rate = stat.mean(rev_rates)
  rev_rate_quants = stat.quantiles(rev_rates, .25, .75)
  
  print "* %d/%d trajectories went back." % (num_rev, num_trials)
  if num_rev > 0:
    print "  Estimated Reverse Collision Rate: %f /uM/s" % (rev_coll_rate / 1000000)
    print "  Estimated Reverse Trajectory Rate: %f (50%% range: %r) /s" % (rev_rate, rev_rate_quants)

  
if __name__ == "__main__":
  import re
  from optparse import OptionParser, OptionGroup
  
  # Parse command line options.
  usage = "usage: %prog [options] BASENAME"
  parser = OptionParser(usage=usage)
  parser.set_defaults(run_kin=True, cleanup=False)
  #parser.set_defaults(verbose=True)
  #parser.add_option("-v", "--verbose", action="store_true", dest="verbose")
  # TODO: implement quiet
  #parser.add_option("-q", "--quiet", action="store_false", dest="verbose")
  parser.add_option("--save", help="Saved state file [defaults to BASENAME.save]", metavar="FILE")
  parser.add_option("--design", help="Design file [defaults to BASENAME.mfe]", metavar="FILE")
  parser.add_option("--seqs", help="Sequences output file [defaults to BASENAME.seqs]", metavar="FILE")
  parser.add_option("--strands", help="Produce a strands-to-order file", metavar="FILE")
  #TODO: parser.add_option("--kinetic", help="Custom kinetics output file, defaults to BASENAME.kin")
  
  kin_parser = OptionGroup(parser, "Kinetics Options")
  kin_parser.add_option("--no-kin", action="store_false", dest="run_kin", help="Don't run kinetics")
  kin_parser.add_option("--trials", type="int", default=24, help="Number of trials to run [Default = %default]")
  kin_parser.add_option("--time", type="float", default=100000, help="Simulation seconds [Default = %default]")
  kin_parser.add_option("--temp", type="float", default=25.0, help="Degrees Celcius [Default = %default]")
  kin_parser.add_option("--conc", type="float", default=1.0, help="Concentration for all molecules (uM) [Default = %default]")
  
  kin_parser.add_option("--spurious", action="store_true", help="Run pairwise kinetic simulations to look for spurious interaction.")
  kin_parser.add_option("--spurious-time", type="float", default=10.0, help="Simulation seconds for spurious simulation [Default = %default]")
  
  kin_parser.add_option("--no-cleanup", action="store_false", dest="cleanup", help="Keep temporary files. [Default temporarily]")
  kin_parser.add_option("--cleanup", action="store_true", dest="cleanup", help="Remove temporary files after use.")
  parser.add_option_group(kin_parser)
  
  (options, args) = parser.parse_args()
  
  # Get basename
  if len(args) < 1:
   parser.error("missing required argument BASENAME")
  basename = args[0]
  # Infer the basename if a full filename is given
  p = re.match(r"(.*)\.(save|mfe)\Z", basename)
  if p:
    basename = p.group(1)
  
  # Set filename defaults
  if not options.save:
    options.save = basename + ".save"
  if not options.design:
    options.design = basename + ".mfe"
  if not options.seqs:
    options.seqs = basename + ".seqs"
  
  finish(options.save, options.design, options.seqs, options.strands, options.run_kin, options.cleanup, options.trials, options.time, options.temp, options.conc, options.spurious, options.spurious_time)
