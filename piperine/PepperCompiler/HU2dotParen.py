from sys import version_info as v
assert (2,3) <= v < (3,), "Python >= 2.3 required (python 3 not tested)"

import string

import HU_grammar as hu_gram
import DotParen_grammar as dp_gram
import exDotParen_grammar as xdp_gram

def HU2dotParen(hu):
  """Convert Zadeh's HU notation to dot-paren notation."""
  ## private recursive function expand
  def expand(p):
    dotParen = ""
    for term in p:
      if term[0] == "+":
        dotParen += "+"
      elif term[0] == "U":
        dotParen += "."*term[1]
      elif term[0] == "H":
        dotParen += "("*term[1]
        dotParen += expand(term[2])
        dotParen += ")"*term[1]
    return dotParen
  ## end of nested function
  par = hu_gram.parse(hu) # parse expression
  return expand(par)   # expand it in dotParen notation

def extended2dotParen(sin):
  """Convert extended dot-paren to standard dot-paren notation."""
  # Parse the extended notation
  p = xdp_gram.parse(sin)
  
  sout = []
  # For each term ...
  for num, symb in p:
    # ... put 'num' repetitions of 'symb'
    sout += symb * num
  sout = string.join(sout, "")
  try:
    dp_gram.parse(sout)
  except xdp_gram.ParseException:
    raise Exception, "Parens don't match %s -> %s" % (sin, sout)
  return sout

def count_dots(p):
  num = 0
  try:
    term = p.next()
    while term == ".":
      num += 1
      term = p.next()
    return num, term
  except StopIteration:
    return num, None

def count_parens(expr):
  num = 0
  while len(expr) == 1 and type(expr[0]) is not str:
    num += 1
    expr = expr[0][1]
  return num, expr

def dotParen2HU(dotParen):
  ## private recursive function to resolve parsed text
  def resolve(p):
    p = iter(p)
    hu = ""
    try:
      term = p.next()
      while term:
        if term == "+":
          hu += "+ "
          term = p.next()
        elif term == ".":
          num, term = count_dots(p)
          hu += "U%d " % (num + 1)
        else:
          assert term[0] == "("
          num, expr = count_parens(term[1])
          hu += "H%d(" % (num + 1)
          hu += resolve(expr)
          hu += ") "
          term = p.next()
    except StopIteration:
      pass
    return hu.strip()
  ## end of nested function
  par = dp_gram.parse(dotParen)
  return resolve(par)
