"""Selector for DNAfold."""
import sys

from utils import error

try:
  import config_choices as config
except ImportError:
  error("DNA Circuit Compiler is not configured, please run config.py")

def DNAfold(seq, temp=25):
  """Run the installed thermodynamic mfe package"""
  if config.thermo == config.NUPACK:
    import DNAfold_Nupack
    return DNAfold_Nupack.DNAfold(seq, temp, exe=config.nupack_path)
  else:
    import DNAfold_Vienna
    return DNAfold_Vienna.DNAfold(seq, temp, exe=config.vienna_path, par_file=config.par_file)
