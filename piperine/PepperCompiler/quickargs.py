"""
Quickly grabs arguments from the command line and uses them as inputs (possibly
keyworded of a function).
"""
import sys
import string
import re

def arg_eval(arg):
  """Evaluate an argument."""
  try:
    return eval(arg)
  except: # If it doesn't evaluate to a real python type assume it's a string.
    return arg

def key_eval(arg):
  """Split a keyword command line argument."""
  assert arg.count("=") == 1, 'Argument "%s" is invalid' % arg
  name, val = string.split(arg, "=")
  return name, arg_eval(val)

def get_args(argv):
  """Read arguments with python syntax, evaluating each string before sending it."""
  # Keywords are of the form name=value
  keys = dict([key_eval(arg) for arg in argv if "=" in arg])
  # Arguments are just values
  args = [arg_eval(arg) for arg in argv if "=" not in arg]
  
  return args, keys

bad_call = re.compile(r"(takes at (least)|(most) \d* arguments)|(got an unexpected keyword argument)")
def call(f):
  """Call f with args supplied by command line."""
  args, keys = get_args(sys.argv[1:])
  
  try:
    return f(*args, **keys)
  except TypeError, e:
    if re.search(bad_call, e.message):
      print e.message
      print "Incorrect usage."
      help(f)
    else:
      raise


def test_printargs(*args, **keys):
  print "Arguments are:", args
  print "Keys are:", keys

def test_1arg(x):
  """
  test_1arg(x) = x
  
  Must be called with exactly one argument "x".
  """
  return x

if __name__ == "__main__":
  call(test_printargs)
  call(test_1arg)
