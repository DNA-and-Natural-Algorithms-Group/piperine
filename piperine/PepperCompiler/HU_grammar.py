"""Grammar for Joe Zadeh's Helix Unpaired notation."""

import sys
sys.path += ("..",) # Extend path to find pyparsing.py
from pyparsing import *

Map = lambda func: (lambda s, t, l: map(func, l))

int_ = Word(nums).setParseAction(Map(int))

# This is a recursive grammar and requires using the 'exression' before it is defined
expr = Forward()

strand_break = Group("+")
unpaired = Group("U" - int_)
helix = Group("H" - int_ + Suppress("(") + Group(expr) + Suppress(")"))

term = strand_break | unpaired | helix

# Now we can define expression which is recursively defined
expr << ZeroOrMore(term)

def parse(s):
  """Parse a string in extended dot-paren notation."""
  return expr.parseString(s, parseAll=True)

# For Example: U6 H7(U4 +) U3 -> [["U", 6], ["H", 7, [["U", 4], ["+"]]], ["U", 3]]
