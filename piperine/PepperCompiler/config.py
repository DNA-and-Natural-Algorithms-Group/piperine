#!/usr/bin/env python
"""
Configure the DNACircuitCompiler to run on your system.
"""

# TODO: configure more things.
#  * spuriousC
#  * Multistrand
#  * ...

import sys
import os

from utils import search_file, search_sys_path, red, green

NUPACK = "NUPACK"
VIENNA = "VIENNA"

# TODO: purhaps, these could be more intelligent
is_nupack_mfe = is_vienna_rnafold = is_vienna_par_file = os.path.isfile

def find_nupack_mfe():
  if "NUPACKHOME" in os.environ:
    path = os.path.join(os.environ["NUPACKHOME"], "bin", "mfe")
    if is_nupack_mfe(path):
      return path
    else:
      print "Your NUPACKHOME environment variable is incorrect (%s)." % os.environ["NUPACKHOME"]
  
  # Search for executable in system path
  return search_sys_path("mfe")

def find_vienna_rnafold():
  if "VIENNAHOME" in os.environ:
    path = os.path.join(os.environ["VIENNAHOME"], "bin", "RNAfold")
    if is_vienna_rnafold(path):
      return path
    path = os.path.join(os.environ["VIENNAHOME"], "Progs", "RNAfold")
    if is_vienna_rnafold(path):
      return path
  
  # Search for executable in system path
  return search_sys_path("RNAfold")

def find_vienna_par_file():
  # If users have set the nice $VIENNAHOME environment variable, it's easy.
  if "VIENNAHOME" in os.environ:
    path = os.path.join(os.environ["VIENNAHOME"], "dna.par")
    if is_vienna_par_file(path):
      return path
  
  # Otherwise search the system path
  else:
    dirs = os.environ["PATH"].split(os.path.pathsep)
    dirs += [".", "~"]
    return search_file("dna.par", dirs)


print "This configuration script will set up the DNA Circuit Compiler to run on your system."

print "Please choose which thermodynamics package to use:"
print "(1) Nupack mfe (recomended/default)"
print "(2) Vienna RNAfold"
thermo = raw_input("?")

# Nupack mfe
if thermo in ("", "1", "(1)"):
  thermo = NUPACK
  print "Searching for Nupack mfe"
  nupack_path = find_nupack_mfe()
  if nupack_path:
    print green("Found Nupack mfe") + " at %s (version unknown)" % nupack_path
    other = raw_input("To use a different version, enter the path [enter nothing for default above]:")
    if other:
      # TODO: Don't exit
      assert is_nupack_mfe(other), "Bad path"
      nupack_path = other
  else:
    print red("Could not find Nupack mfe")
    nupack_path = raw_input("Please enter path to Nupack mfe:")
    # TODO: Don't exit
    assert is_nupack_mfe(nupack_path), "Bad path"

# Vienna RNAfold
elif thermo in ("2", "(2)"):
  thermo = VIENNA
  print "Searching for Vienna RNAfold"
  vienna_path = find_vienna_rnafold()
  if vienna_path:
    print green("Found Vienna RNAfold") + " at %s (version unknown)" % vienna_path
    other = raw_input("To use a different version, enter the path [enter nothing for default above]:")
    if other:
      # TODO: Don't exit
      assert is_vienna_rnafold(other), "Bad path"
      vienna_path = other
  else:
    print red("Could not find Vienna RNAfold")
    vienna_path = raw_input("Please enter path to Vienna RNAfold:")
    # TODO: Don't exit
    assert is_vienna_rnafold(vienna_path), "Bad path"
  
  print "Searching for dna.par file"
  par_file = find_vienna_par_file()
  if par_file:
    print green("Found Vienna parameter file") + " at %s (version unknown)" % par_file
    other = raw_input("To use a different version, enter the path [enter nothing for default above]:")
    if other:
      # TODO: Don't exit
      assert is_vienna_par_file(other), "Bad path"
      par_file = other
  else:
    print red("Could not find Vienna parameter file")
    par_file = raw_input("Please enter path to Vienna parameter file:")
    # TODO: Don't exit
    assert is_vienna_par_file(par_file), "Bad path"

else:
  # TODO: Don't exit
  raise Exception, "Unexpected choice"




## Write configuration file
# We need to write the configuration file in the same directory as this script.
path = sys.argv[0]  # path to this script
dir = os.path.dirname(path)  # directory this script is in
config_file = os.path.join(dir, "config_choices.py")

f = open(config_file, "w")
f.write("NUPACK = 'NUPACK'\n")
f.write("VIENNA = 'VIENNA'\n")
if thermo == NUPACK:
  f.write("thermo = NUPACK\n")
  f.write("nupack_path = %r\n" % nupack_path)
else:
  f.write("thermo = VIENNA\n")
  f.write("vienna_path = %r\n" % vienna_path)
  f.write("par_file = %r\n" % par_file)
f.close()
