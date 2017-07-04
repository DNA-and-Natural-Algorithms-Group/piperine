#!/usr/bin/env python

import sys
import itertools as it

num_left  = int(sys.argv[1])
num_right = int(sys.argv[2])

left_indices  = range(1, num_left + 1)
right_indices = range(1, num_right + 1)

filename = 'Seesaw_%d_%d.comp' % (num_left, num_right)

f = open(filename, 'w')

# Declaration statement
f.write('declare component Seesaw_%d_%d(toe, rec, clamp): ' % (num_left, num_right))
f.write(' + '.join(['left%d' % (i) for i in left_indices]))
f.write(' -> ')
f.write(' + '.join(['right%d' % (i) for i in right_indices]))

# Domains
f.write('\n\n### Domains\n')
f.write('sequence toe  = "?N" : <toe>   # Universal toehold domain\n')
f.write('sequence base = "?N" : <rec>   # Base recognition domain on seesaw gate\n')
f.write('sequence c    = "?N" : <clamp> # Clamp domain\n')
for i in left_indices:
  f.write('sequence l%d   = "?N" : <rec>   # Left recognition domain\n' % (i))
for i in right_indices:
  f.write('sequence r%d   = "?N" : <rec>   # Right recognition domain\n' % (i))

# I/O connectors
f.write('\n### I/O connectors\n')
for i in left_indices:
  f.write('sequence left%d  = c l%d c toe c base c : <4*clamp + 2*rec + toe>\n' % (i, i))
for i in right_indices:
  f.write('sequence right%d = c base c toe c r%d c : <4*clamp + 2*rec + toe>\n' % (i, i))

# Strands
f.write('\n### Strands\n')
f.write('strand Base   = c* toe* c* base* c* toe* c* : <4*clamp + 2*toe + rec>  # Base of the seesaw itself\n')
for i in left_indices:
  f.write('strand Left%d  = c l%d c toe c base c : <4*clamp + 2*rec + toe>  # Strand which fits on left side of seesaw\n' % (i, i))
for i in right_indices:
  f.write('strand Right%d = c base c toe c r%d c : <4*clamp + 2*rec + toe>  # Strand which fits on right side of seesaw\n' % (i, i))

# Structures
f.write('\n### Structures\n')
f.write('# Single-stranded:\n')
for i in left_indices:
  f.write('structure LEFT%d  = Left%d  : domain .......  # A single-stranded left strand\n' % (i, i))
for i in right_indices:
  f.write('structure RIGHT%d = Right%d : domain .......  # A single-stranded right strand\n' % (i, i))
f.write('\n# Double-stranded:\n')
for i in left_indices:
  f.write('structure SEESAW_LEFT%d  = Base + Left%d  : domain ..(((((+..)))))  # Seesaw in the tipped left position\n' % (i, i))
for i in right_indices:
  f.write('structure SEESAW_RIGHT%d = Base + Right%d : domain (((((..+)))))..  # Seesaw in the tipped right position\n' % (i, i))

# Reactions
f.write('\n### Reactions\n')
for i, j in it.product(left_indices, right_indices):
  f.write('kinetic LEFT%d + SEESAW_RIGHT%d -> SEESAW_LEFT%d + RIGHT%d  # Forward reaction\n' % (i, j, i, j))
  f.write('kinetic SEESAW_LEFT%d + RIGHT%d -> LEFT%d + SEESAW_RIGHT%d  # Reverse reaction\n' % (i, j, i, j))
  f.write('\n')

