import re

import system_class
import component_class

def split(string):
  #return string.split()
  return re.findall(r"[\w?*]+", string)

def parse_declare(line):
  # TODO: allow comments and variable spacing. Maybe use re.VERBOSE
  p = re.match(r"""declare\ (system|component)\ 
                   (\w+)              # Declare name
                   (?:\((.*?)\))?:    # Parameter names
                   (.*) -> (.*)\n     # Ins and outs""", line, re.VERBOSE)
  assert p, "First line of file must declare system/component.\n" + line # TODO: add syntax
  type_, name, param_names, ins, outs = p.groups("")
   
  # Expand out the lists
  param_names = split(param_names)
  ins  = split(ins)
  outs = split(outs)
  
  return type_, name, param_names, ins, outs

# TODO: Make into a file-like object for efficiency and preserving original file.
def preprocess(f, params):
  """Preprocesses a file. Applies parameters and removes comments and newlines."""
  doc = []
  for line in f:
    # Strip comments
    line = re.sub(r"#.*\n", r"\n", line)
    
    # Skip empty lines
    if re.match(r"\s*\Z", line):
      continue
    
    # Replace all <expression> with eval(expression) (Caution: security hole)
    def eval_brackets(s):
      return str(eval(s.group(1), params))
    line = re.sub(r"<(.*?)>", eval_brackets, line)
    
    doc.append(line)
  return doc

def load_file(filename, args):
  """Load a system or component from filename."""
  f = open(filename, "r")
  line = f.readline()
  
  # Read first line.
  type_, name, param_names, ins, outs = parse_declare(line)
  # TODO: check that name matches filename
  assert len(args) == len(param_names), "Argument length mismatch.\n%d Params: %r\n%d Args: %r\n" % (len(param_names), param_names, len(args), args) # TODO: more info
  
  # Do parameter substitution.
  params = dict(zip(param_names, args))
  doc = preprocess(f, params)
  
  print doc
  
  if type_ == "system":
    return load_system(doc, name, ins, outs)
  else: # type_ == "component"
    return load_component(doc, name, ins, outs)


## System stuff
def parse_import(rest):
  # TODO: allow variable spacing. Maybe use re.VERBOSE
  p = re.match(r"(\w+)(?: as (\w+))?\n", rest)
  assert p, "Import syntax incorrect.\n" + rest # TODO: add line and syntax
  path, name = p.groups()
  
  return path, name

def parse_component(rest):
  # TODO: allow variable spacing. Maybe use re.VERBOSE
  p = re.match(r"(\w+) = (\w+)(?:\((.*?)\))?: (.*) ?-> ?(.*)\n", rest)
  assert p, "Component syntax incorrect" # TODO: add line and syntax
  component_name, templ_name, templ_args, ins, outs = p.groups("")
  
  # Expand out the lists
  templ_args = split(templ_args)
  ins  = split(ins)
  outs = split(outs)
  
  return component_name, templ_name, templ_args, ins, outs

def load_system(doc, name, ins, outs):
  """Build a system object from the commands in the file."""
  system = system_class.System(name, None, ins, outs)
  for line in doc:
    command, rest = line.split(None, 1) # Split off first word # TODO: Put in try/except for when this fails
    if command == "import":
      path, name = parse_import(rest)
      system.add_import((path, name)) # TODO: Check rest ... and fix
    elif command == "component":
      component_name, templ_name, templ_args, ins, outs = parse_component(rest)
      system.add_component(component_name, templ_name, templ_args, ins, outs)
    else:
      assert False, "Command " + command + " not defined in system." # TODO: where
  return system


## Component stuff
def parse_seq(rest):
  # TODO: allow variable spacing. Maybe use re.VERBOSE
  p = re.match(r"(\w+) = (.+) : (\d+)\n", rest)
  assert p, "Sequence syntax incorrect.\n" + `rest` # TODO: add line and syntax
  name, const, length = p.groups("")
  
  length = int(length)
  #Process constraint list
  const_temp = split(const)
  const = []
  for item in const_temp:
    p = re.match(r"(\d+|\?)?(\w)", item) # TODO: Catch parse error
    num, code = p.groups(1) # No number means 1
    const.append((num, code))
  
  return name, const, length

def parse_strand(rest):
  # TODO: allow variable spacing. Maybe use re.VERBOSE
  p = re.match(r"(\w+) = (.+) : (\d+)\n", rest)
  assert p, "Strand syntax incorrect.\n" + `rest` # TODO: add line and syntax
  name, const, length = p.groups("")
  
  length = int(length)
  #Process constraint list
  const_temp = split(const)    # TODO: More processing necessary
  const = []
  for item in const_temp:
    p = re.match(r"(\d+|\?)?([ACGTUNSWRYMKVHBD])", item)
    if p: # It's an anonymous constraint
      num, code = p.groups(1) # No number means 1
      const.append((num, code))
    else: # It's a sequence name
      p = re.match(r"(\w+)(\*)?", item)
      seq_name, wc = p.groups("")
      const.append(("Sequence", (seq_name, wc)))
  
  return name, const, length

def load_component(doc, name, ins, outs):
  """Build a component object from the commands in the file."""
  # TODO: fix
  ins = [(x, None) for x in ins]
  outs = [(x, None) for x in outs]
  component = component_class.Component(name, None, ins, outs)
  for line in doc:
    command, rest = line.split(None, 1) # Split off first word
    if command == "sequence":
      name, const, length = parse_seq(rest)
      component.add_sequence(name, const, length)
    # TODO: Add super sequence
    elif command == "strand":
      name, const, length = parse_strand(rest)
      component.add_strand(name, const, length)
    # TODO: rest of commands
    else:
      assert False, "Command " + command + " not defined in component.\n" + line # TODO: where
  return component

if __name__ == "__main__":
  import sys
  
  component_class.DEBUG = system_class.DEBUG = True
  
  filename = sys.argv[1]
  args = map(eval, sys.argv[2:])
  print load_file(filename, args)
