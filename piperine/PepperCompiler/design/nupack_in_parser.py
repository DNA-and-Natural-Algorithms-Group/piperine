"""Nucleic Acid template design grammar"""

from nupack_in_class import Spec

from ..pyparsing import *

## Some globals
# Pyparsing shortcuts
K = CaselessKeyword
S = Suppress
O = Optional
H = Hidden = lambda x: Empty().setParseAction(lambda s,t,l: x)  # A hidden field, it tags an entry
List = lambda x: Group(OneOrMore(x))
Map = lambda func: (lambda s,l,t: map(func, t) )
# syntax names and sets
# Codes for different subsets of nucleotides
NAcodes = "ACGTUNSWRYMKVHBD"
struct = "structure"
seq = "sequence"
prevent = "prevent"
app = "apply"
obj = "objective"

ParserElement.setDefaultWhitespaceChars(" \t")


## Define Grammar
var = Word(alphas+"_-", alphanums+"_-") # Variable name
var_list = List(var)

integer = Word(nums).setParseAction(Map(int))
float_ = Word(nums+"-.eE").setParseAction(Map(float))

sequence = Word(NAcodes)

seq_var = Group(var + Optional("*", default=""))
seq_list = List(seq_var)

### TODO-maybe: Actually parse it.
secondary_struct = Word( nums+".()+ " )

# structure <name> = <secondary structure>
struct_stat = K(struct) + var + S("=") + secondary_struct
# sequence <name> = <constraints>
seq_stat  = K(seq)  + var + S("=") + sequence
# prevent <list of structs> < <float>
prevent_stat = K(prevent) + var_list + S("<") + float_
# <struct name> : <list of seqs>
apply_stat = H(app) + var + S(":") + seq_list
# <struct name> < <float>
obj_stat = H(obj) + var + S("<") + float_

statement = struct_stat | seq_stat | prevent_stat | apply_stat | obj_stat

document = StringStart() + delimitedList(O(Group(statement)), "\n") + StringEnd()
document.ignore(pythonStyleComment)

def load_design(filename):
  """Load component design file"""
  try:
    # Load data
    statements = document.parseFile(filename)
  
  except ParseException, e:
    print
    print "Parsing error in template:", filename
    print e
    sys.exit(1)
  
  # Build data
  spec = Spec()
  for stat in statements:
    #print list(stat)
    if stat[0] == struct:
      spec.add_structure(*stat[1:])
    elif stat[0] == seq:
      spec.add_sequence(*stat[1:])
    elif stat[0] == app:
      spec.add_apply(*stat[1:])
  return spec
