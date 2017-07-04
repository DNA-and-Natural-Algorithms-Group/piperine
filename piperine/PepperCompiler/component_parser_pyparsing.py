"""DNA component design grammar with pyparsing"""

import sys

from HU2dotParen import extended2dotParen, HU2dotParen
from var_substitute import process
from utils import print_linenums, error

from pyparsing import *

## Some globals
# Pyparsing shortcuts
K = CaselessKeyword
S = Suppress
O = Optional
H = Hidden = lambda x: Empty().setParseAction(lambda s, l, t: x)  # A hidden field, it tags an entry
Map = lambda func: (lambda s, l, t: map(func, t) )

def List(expr, delim=""):
  """My delimited list. Allows for length zero list and uses no delimiter by default."""
  if not delim:
    return Group(ZeroOrMore(expr))
  else:
    return Group(Optional(expr + ZeroOrMore( Suppress(delim) + expr)))

def Flag(expr):
  """A flag identifier. It is either present or not and returns True or False."""
  p = Optional(expr)
  p.setParseAction(lambda s, l, t: bool(t))
  return p

import string
lowers = string.lowercase


nucleotide_flag = "Anonymous"
sequence_flag = "Sequence"
domains_flag = "domains"
def result2list(foo):
  """Convert from ParseResults to normal list."""
  if isinstance(foo, ParseResults):
    return [result2list(bar) for bar in foo]
  else:
    return foo


# Syntax names and sets
# Codes for different subsets of nucleotides
NAcodes = "ACGTUNSWRYMKVHBD"

decl = "declare"
component = "component"
seq = "sequence"
sup_seq = "super-sequence"
strand = "strand"
struct = "structure"
kin = "kinetic"

# Don't ignore newlines!
ParserElement.setDefaultWhitespaceChars(" \t")


## Define Grammar
var = Word(alphas, alphanums+"_") # Variable name
integer = Word(nums).setParseAction(Map(int))
float_ = Word(nums+"+-.eE").setParseAction(Map(float))

# Sequence const could be ?N, 3N or N
seq_const = S('"') + List(Group( Optional(integer | "?", default=1) +
                                 Word(NAcodes, exact=1))) + S('"')

seq_name = var
seq_var = Group(seq_name + Flag("*"))

# Strand definition could be:
#  1) Some basic sequence constraints like 5N or 2S or
#  2) A sequence variable (possibly complemented with *) (Must start with lowercase letter)
sup_seq_const = ( Group(H(nucleotide_flag) + seq_const) |
                  Group(domains_flag + S("(") + seq_var + S(")")) |
                  Group(H(sequence_flag) + seq_var) )

strand_var = var
struct_var = var

# Signals used in the declare line seq
signal_var = Group(seq_var + Optional(S("(") + struct_var + S(")"), default=None))

# Secondary structure can be dot-paren, extended dot-paren or HU notation.
# 'domain' keyword means that the helix/unpaired segments are by domain
## TODO: deal with errors. Super-confusing error messages right now.
exDotParen = Word(nums + ".()+ " ).setParseAction(Map(extended2dotParen))
HUnotation = Word(nums + "UH()+ ").setParseAction(Map(HU2dotParen))
secondary_struct = Group( Flag("domain") + (exDotParen | HUnotation) )


# declare component <component name>(<params>): <inputs> -> <outputs>
params = O(S("(") + List(var, ",") + S(")"), default=[])
decl_stat = K(decl) + S(component) + var + params + S(":") + List(signal_var, "+") + S("->") + List(signal_var, "+")
def parse_declare_statement(statement):
  return result2list(decl_stat.parseString(statement, parseAll=True))[1:]

# sequence <name> = <constraints> : <length>
seq_stat = K(seq)  + seq_name + S("=") + List(sup_seq_const) + \
           Optional(S(":") + integer, default=None)
def parse_general_sequence_statement(statement):
  return result2list(seq_stat.parseString(statement, parseAll=True))[1:]

# super-sequence <name> = <constraints> : <length>
sup_seq_stat = K(sup_seq)  + seq_name + S("=") + List(sup_seq_const) + \
           Optional(S(":") + integer, default=None)
def parse_super_sequence_statement(statement):
  return result2list(seq_stat.parseString(statement, parseAll=True))[1:]

# strand <name> = <constraints / sequences> : <length>
strand_stat = K(strand) + Flag("[dummy]") + strand_var + S("=") + List(sup_seq_const) + \
              Optional(S(":") + integer, default=None)
def parse_strand_statement(statement):
  return result2list(strand_stat.parseString(statement, parseAll=True))[1:]

# structure <optinoal opt param> <name> = <strands> : <secondary structure>
opt = Optional(   K("[no-opt]").setParseAction(lambda s, t, l: False) | \
                ( Suppress("[") + float_ + Suppress("nt]") ),
                  default=1.0)
struct_stat = K(struct) + opt + struct_var + S("=") + List(strand_var, "+") + S(":") + secondary_struct
def parse_structure_statement(statement):
  return result2list(struct_stat.parseString(statement, parseAll=True))[1:]

# kin <inputs> -> <outputs>
kin_stat = K(kin) + List(struct_var, "+") + S("->") + List(struct_var, "+")
def parse_kinetic_statement(statement):
  return result2list(kin_stat.parseString(statement, parseAll=True))[1:]


statement = seq_stat | strand_stat | struct_stat | kin_stat
document = StringStart() + ZeroOrMore(S("\n")) + Group(decl_stat) + S("\n") + \
           List(O(Group(statement)), "\n") + StringEnd()
document.ignore(pythonStyleComment)



def load_component(filename, args, prefix):
  """Load component file"""
  from component_class import Component
  try:
    # Open file and do parameter substitution
    doc = substitute(filename, args)
  except ParseBaseException, e:
    print
    print "Parsing error in component:", filename
    print e
    sys.exit(1)
    
  try:
    # Load data
    declare, statements = result2list(document.parseString(doc, parseAll=True))
  except ParseBaseException, e:
    print
    print_linenums(doc)
    print "Parsing error in component:", filename
    print e
    sys.exit(1)
  
  command, name, params, inputs, outputs = declare
  # Build data
  component = Component(name, prefix, params)
  for statement in statements:
    #print list(stat)
    if statement[0] == seq:
      command, name, const, length = statement
      if len(const) == 1 and const[0][0] == "Anonymous":
        component.add_sequence(name, const[0][1], length)
      else:
        component.add_super_sequence(name, const, length)
    
    elif statement[0] == strand:
      component.add_strand(*statement[1:])
    
    elif statement[0] == struct:
      component.add_structure(*statement[1:])
    
    elif statement[0] == kin:
      component.add_kinetic(*statement[1:])
    
    else:
      raise Exception, "Unexpected statement:\n%s" % statement
  component.add_IO(inputs, outputs)
  return component

def substitute(filename, args):
  # Parse for function declaration
  param_names = decl_stat.parseFile(filename)[2]
  params = {}
  if len(param_names) != len(args):
    error("Component %s takes %d parameters %r, %d given %r." % (filename, len(param_names), tuple(param_names), len(args), tuple(args)))
  for name, val in zip(param_names, args):
    params[name] = val
  return process(filename, params)
