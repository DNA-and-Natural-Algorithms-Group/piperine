import os

class BadFilename(Exception):
  """File not found"""
  pass

def find_file(name, ext=".des"):
  """Find the actual filename, which is either name or name + ext """
  if os.path.isfile(name):
    return name
  elif os.path.isfile(name + ext):
    return name + ext
  else:
    raise BadFilename, "File not found: neither %s nor %s exist." % (name, name + ext)
