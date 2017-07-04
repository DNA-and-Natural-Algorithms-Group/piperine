"""
Does some simple processing:
 * <expr>  are evaluated
 * {a,b} are separated onto multiple lines (as in shell scripting)
 
Lines begining with length may be used to set variables

So, the string:
 I have {two,three,four} apples worth $<(3-5+17)/4>.
goes to
 I have two apples worth $3.
 I have three apples worth $3.
 I have four apples worth $3.
"""

import string
import re

def process_filename(infilename, params):
  return process_list(open(infilename, "r"), params)

process = process_filename

def process_list(lines, params):
  if "__builtins__" not in params:
    params["__builtins__"] = None
  
  out = ""
  for line in lines:
    # Remove comments
    line = re.sub(r"#.*", "", line)
    
    # Evaluate length lines
    m = re.match(r"\s*length\s+(\w+)\s*=\s*(.*)", line)
    if m:
      name, val = m.groups()
      val = eval(val, params)
      params[name] = val
      continue
    
    def eval_brackets(s):
      return str(eval(s.group(1), params))
    # Replace all <expression> with eval(expression) (Caution: security hole)
    line = re.sub(r"<([^<>]*?)>", eval_brackets, line)
    
    def duplicate(line):
      # Replace first {...,expr,...} with expr for each expression in list
      m = re.search(r"{([^{}]*?)}", line)
      if m:
        ops = m.group(1).split(",")
        start = line[:m.start()]
        end = line[m.end():]
        this_out = ""
        for op in ops:
          this_out += duplicate(start+op+end)
        return this_out
      else:
        return line
    lines = duplicate(line)
    
    # Don't add blank lines
    if lines.strip():
      out += lines
  
  return out


if __name__ == "__main__":
  import sys
  
  import system_parser
  import component_parser
  
  def substitute(filename, args):
    # Find first statement
    for line in open(filename, "r"):
      line = re.sub(r"#.*\n", "", line)
      line = line.strip()
      if line:
        break
    # Parse as system/component declaration
    if re.match(r".*\.sys\Z", filename):
      param_names = system_parser.decl_stat.parseFile(filename)[2]
    elif re.match(r".*\.comp\Z", filename):
      param_names = component_parser.parse_declare_statement(line)[1]
    else:
      raise ValueError, "File %s is neither system nor component type." % filename
    # and do substitution.
    params = {}
    assert len(param_names) == len(args), (param_names, args)
    for name, val in zip(param_names, args):
      params[name] = val
    return process_filename(filename, params)

  filename = sys.argv[1]
  args = map(eval, sys.argv[2:])
  print substitute(filename, args)
