"""Useful classes and functions."""
from __future__ import division

import math
import os
import random
import re
import string
import sys
import tempfile

DEBUG = False  # For debugging and unittests, set to true to raise exception on error and warning instead of sys.exit(1) and nothing.

def print_linenums(text):
   """Print text with line numbers prepended. Start counting at line 1 to fit with pyparsings line numberings"""
   for n, line in enumerate(text.split("\n")):
     line_num = "%3d:" % (n + 1)
     print line_num + line

## Custom regular expression wrapper
class ParseException(Exception): pass

def match(regex, line):
  """Match *entire* line to regex converting all spaces into '\s+' and allowing trailing spaces."""
  regex = regex.replace(" ", r"\s+")
  regex += r"\s*\Z"
  return re.match(regex, line)

## Search for file in list of directories
def search_file(filename, search_path):
  """Find a file with this name in search_path."""
  for dir in search_path:
    file_path = os.path.join(dir, filename)
    if os.path.isfile(file_path):
      return file_path
  return None

def search_sys_path(filename):
  """Search for file in system path."""
  search_path = os.environ["PATH"].split(os.path.pathsep)
  return search_file(filename, search_path)


## ANSI Colors
red_tag = "\033[31;1m"   # ANSI color code for bold red
green_tag = "\033[32;1m" # ANSI color code for bold green
reset_tag = "\033[m"     # ANSI color code for reset color

def red(text):
  return red_tag + text + reset_tag
def green(text):
  return green_tag + text + reset_tag


## Error and warning messages
def error(text):
  """Print a formatted error message and exit."""
  if not DEBUG:
    sys.stderr.write(red("ERROR: %s\n" % text))
    sys.exit(1)
  else:
    raise Exception, red("ERROR: %s\n" % text)

def warning(text):
  """Print a formatted warning message."""
  if not DEBUG:
    sys.stderr.write(red("Warning: %s\n" % text))
  else:
    raise Exception, red("Warning: %s\n" % text)


def mktemp(mode, *args, **keys):
  """Creates a temporary file. Returns the file and filename.
     Like tempfile.mkstemp except that it actually creates a file object."""
  fd, filename = tempfile.mkstemp(*args, **keys)
  file_ = os.fdopen(fd, mode)
  return file_, filename


def _dummy(*args, **kw):
  """A method not available function."""
  raise Exception, "methods not available"

## Generic Objects
class ordered_dict(dict):
  """A standard dictionary that remembers the order you added items in.
     Only supports __setitem__, __iter__, keys, values, items and read-only ops."""
  def __init__(self):
    self.order = []
    dict.__init__(self)
  
  def get_index(self, index):
    return self[self.order[index]]
  
  def __setitem__(self, key, value):
    if key not in self:
      self.order.append(key)
    dict.__setitem__(self, key, value)
  def __iter__(self):
    for key in self.order:
      yield key
  def keys(self):
    return self.order
  def values(self):
    return [self[key] for key in self]
  def items(self):
    return [(key, self[key]) for key in self]
  
  def __delitem__(self, key):
    self.order.remove(key)
    dict.__delitem__(self, key)
  ### TODO-maybe: impliment clear, copy, iter*, ...
  clear = copy = iteritems = iterkeys = itervalues = pop \
        = popitem = update = fromkeys = _dummy

class default_ordered_dict(ordered_dict):
  """An ordered dictionary automatically defaults to a given value on unset key.
     With the call=True flag, the default will be assumed a function and called each time."""
  def __init__(self, default=None, call=False):
    self.default = default
    self.call = call
    ordered_dict.__init__(self)
  def __getitem__(self, key):
    """Get key. If it doesn't exist, set it to default first."""
    if key not in self:
      if not self.call:
        self[key] = self.default
      else:
        self[key] = self.default()
    
    return self.get(key)

class ordered_set(set):
  """A standard set that remembers the order you added items in.
     Only supports add, update, __iter__ and read-only ops."""
  def __init__(self):
    self.order = []
    set.__init__(self)
  def add(self, value):
    if value not in self:
      self.order.append(value)
      set.add(self, value)
  def update(self, other):
    for value in other:
      self.add(value)
  def __iter__(self):
    for value in self.order:
      yield value
  
  ### TODO-maybe: impliment clear, copy, ...
  clear = copy = difference = difference_update = discard = intersection \
        = intersection_update = pop = remove = symmetric_difference \
        = symmetric_difference_update = union = _dummy

class PrintObject(object):
  """Generic default-printable object."""
  def __str__(self):
    attribs = ["%s=%r" % (name, value) for (name, value) in self.__dict__.items()]
    attribs = string.join(attribs, ", ")
    return "%s(%s)" % (self.__class__.__name__, attribs)
  __repr__ = __str__
