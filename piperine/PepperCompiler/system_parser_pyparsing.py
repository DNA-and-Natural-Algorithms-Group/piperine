import sys

from system_class import System
from var_substitute import process
from utils import print_linenums, error

from pyparsing import *

def result2list(foo):
  """Convert from ParseResults to normal list."""
  if isinstance(foo, ParseResults):
    return [result2list(bar) for bar in foo]
  else:
    return foo

## Some globals
K = CaselessKeyword
S = Suppress
O = Optional
Map = lambda func: (lambda s, l, t: map(func, t) )  # A useful mapping function

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

decl = "declare"
system = "system"
import_ = "import"
component = "component"

# Don't ignore newlines!
ParserElement.setDefaultWhitespaceChars(" \t")


## Define Grammar
var = Word(alphas, alphanums+"_") # Variable name
signal = Group(var + Flag("*"))
signal_list = List(signal, "+")

path = Word(alphanums+".-_/~") # Path name in a directory structure
py_chars = printables.replace(",", "").replace(")", "")
python_object = Word(py_chars, py_chars+" ").setParseAction(Map(eval))


# declare system <cicuit name>(<params>): <inputs> -> <outputs>
decl_params = O(S("(") + List(var, ",") + S(")"), default=[])
decl_stat = K(decl) + S(system) + var + decl_params + S(":") + signal_list + S("->") + signal_list
def parse_declare_statement(statement):
  return result2list(decl_stat.parseString(statement, parseAll=True))[1:]

# import Adder, HalfAdder5 as HalfAdder, templates/Crossing_Gates/LastAdder
import_stat = K(import_) + delimitedList(Group(path + O(S("as") + var, default=None)))
def parse_import_statement(statement):
  return result2list(import_stat.parseString(statement, parseAll=True))[1:]

# component <name> = <template name>(<params>): <inputs> -> <outputs>
component_params = O( S("(") + List(python_object, ",") + S(")") , default=[])
component_stat = K(component) + var + S("=") + var + component_params + S(":") + signal_list + S("->") + signal_list
def parse_component_statement(statement):
  return result2list(component_stat.parseString(statement, parseAll=True))[1:]


statement = import_stat | component_stat

document = StringStart() + Group(decl_stat) + S("\n") + \
           List(O(Group(statement)), delim="\n") + StringEnd()
document.ignore(pythonStyleComment)



def load_system(filename, args, prefix, path):
  """Load system connectivity file"""
  try:
    # Open file and do parameter substitution
    doc = substitute(filename, args)
  except ParseBaseException, e:
    print
    print "Parsing error in system:", filename
    print e
    sys.exit(1)
    
  try:
    # Load data
    declare, statements = document.parseString(doc, parseAll=True)
  except ParseBaseException, e:
    print
    print_linenums(doc)
    print "Parsing error in system:", filename
    print e
    sys.exit(1)
  
  x, name, params, inputs, outputs = declare
  # Build data
  system = System(path, name, prefix, params)
  for stat in statements:
    #print list(stat)
    if stat[0] == import_:
      system.add_import(*stat[1:])
    elif stat[0] == component:
      system.add_component(*stat[1:])
    else:
      raise Exception, "Unexpected statement:\n%s" % stat
  system.add_IO(inputs, outputs)
  return system

def substitute(filename, args):
  # Parse for function declaration
  param_names = decl_stat.parseFile(filename)[2]
  params = {}
  if len(param_names) != len(args):
    error("System %s takes %d parameters %r, %d given %r." % (filename, len(param_names), tuple(param_names), len(args), tuple(args)))
  for name, val in zip(param_names, args):
    params[name] = val
  return process(filename, params)
