"""Multistrand wrapper"""
from __future__ import division

import os
import string
import subprocess
import random
import re

from utils import mktemp, error

urandom = random.SystemRandom()
def random_seed():
  """
  Return a random seed x uniform random [0 <= x < 2**31].
  We use urandom to get a more random result than random.getrandbits.
  """
  return hex(urandom.getrandbits(31)).rstrip("L")

# Constants for use in Multistrand
BACK_FLAG = "REVERSE"
FORWARD_FLAG = "FORWARD"

def DNAkinfold(strands, start_struct, back_struct, stop_struct, trials, sim_time, temp, conc, out_interval=-1, cleanup=True):
  """
  strands = dict of strand_name : sequence used in structs
  *_struct = list of structures (each of which is a list of strand_names and a secondary struct)
  temp = temperature (deg C)
  conc = concentration (uM)
  sim_time = max time of sim (approx seconds)
  trials = number of trials to run
  """
  # Note: Assumes that structures are connected complexes
  assert start_struct != stop_struct
  
  f, in_name = mktemp(mode="w", prefix="multi_", suffix=".in")
  out_name = in_name[:-3] + ".out"
  
  # Print input file for Multistrand
  # Strand Definitions
  f.write("#Strands\n")
  for name, seq in strands.items():
    f.write("%s,%s\n" % (name, seq))
  # Start Structure
  f.write("#StartStructure\n")
  for compl in start_struct:
    f.write(string.join(compl.strands, ",") + "\n")  # Seq1,Seq2,Chicken,Seq4
    struct = compl.struct.replace("+", "_")
    f.write(struct + "\n")                        # ((...._..((_))..)..(_..)...)
  f.write("#StopStructures\n")
  # Backward Stop Structure (Where the reaction started and might fall back to)
  for compl in back_struct:
    f.write(string.join(compl.strands, ",") + "\n")  # Seq1,Seq2,Chicken,Seq4
    struct = compl.struct.replace("+", "_")
    f.write(struct + "\n")
  if back_struct: # Don't print the back flag if we have not back structures
    f.write("TAG: %s\n" % BACK_FLAG)
  # Forward Stop Structure (Where the reaction should go)
  for compl in stop_struct:
    f.write(string.join(compl.strands, ",") + "\n")  # Seq1,Seq2,Chicken,Seq4
    struct = compl.struct.replace("+", "_")
    f.write(struct + "\n")
  if stop_struct: # Don't print the forward flag if we have not stop structures
    f.write("TAG: %s\n" % FORWARD_FLAG)
  # Other params
  f.write("#Energymodel=NUPACK_DNA_2_3\n")
  f.write("#Temperature=%f\n" % temp)
  f.write("#Concentration=%f\n" % (conc/1000) ) # Multistrand uses mM instead of uM
  f.write("#SimTime=%f\n" % sim_time)
  f.write("#NumSims=%d\n" % trials)
  f.write("#Logfile=%s\n" % out_name)
  f.write("#OutputInterval=%d\n" % out_interval)   # -1 = Suppress output
  f.write("#StopOptions=2\n")        # Stop on stop structures defined above
  f.write("#SimulationMode=1\n")     # Start with binding two structions and then look at forward and backwards rates.
  # Done
  f.close()
  
  # This'll be a mess if we leave an old outfile there.
  if os.path.isfile(out_name):
    os.remove(out_name)
  # Run Multistrand!
  command = "nice Multistrand < %s" % in_name
  # If we asked for quiet, keep it quiet.
  if out_interval == -1:
    command += " > /dev/null"
  print "$", command
  subprocess.check_call(command, shell=True)
  
  f = open(out_name, "r")
  
  # Note: This will load the entire data of the file into lists (could be memory 
  #       intensive for extremely large data sets).
  forward = []
  reverse = []
  overtime = []
  summary = ""
  for line in f:
    if line[0] == "(":
      parts = line.split()
      assert len(parts) == 5, error("Multistrand format has changed.\n%s" % line)
      coll_rate = float(parts[2])
      flag = parts[3]
      time = float(parts[4])
      if time >= sim_time:
        overtime.append( (time, coll_rate) )
      elif flag == FORWARD_FLAG:
        forward.append( (time, coll_rate) )
      elif flag == BACK_FLAG:
        reverse.append( (time, coll_rate) )
      else:
        assert False, error("Unexpected stop flag '%s' in Multistrand output.\n%s" % (flag, line))
    else:
      summary += line
  f.close()
  
  if cleanup:
  	os.remove(in_name)
  	os.remove(out_name)
  
  return forward, reverse, overtime, summary
