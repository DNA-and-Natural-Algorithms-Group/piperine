"""Grammar for dot-paren notation for secondary structure"""

import sys
sys.path += ("..",) # Extend path to find pyparsing.py
from pyparsing import *

# This is a recursive grammar and requires using the 'exression' before it is defined
expr = Forward()

strand_break = Literal("+")
dot = Literal(".")
parens = Group("(" - Group(expr) + ")")

term = strand_break | dot | parens

# Now we can define expression which is recursively defined
expr << ZeroOrMore(term)

def parse(s):
  """Parse a string in extended dot-paren notation."""
  return expr.parseString(s, parseAll=True)

# For Example: ...((....+))... -> [".", ".", ".", ["(", [["(", [".", ".", ".", ".", "+"], ")"]], ")"], ".", ".", "."]
