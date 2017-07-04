def propagate_constraints(eq, wc):
  """
  Takes a basis of constraints and propagates them to produce all constraints.
  
  Takes dicts eq and wc where 
    y in eq[x] means that x is equal y  and 
    y in wc[x] means that y is complementary to x
  and produces eq_all and wc_all where
    eq_all[x] = {y | y == x} (set of *all* items equal to x)
    wc_all[x] = {y | y ~= x} (set of *all* items equal to x)
  
  NOTE: Assumes that if  y in eq[x]  =>  x in eq[y]  and likewise for wc.
  """
  keys = eq.keys() # Items we are constraining
  assert wc.keys() == keys, "eq.keys() and wc.keys() must be equal to the set of items to constrain over"
  eq_all = {}
  wc_all = {}
  
  for x in keys:
    if x not in eq_all: # If we haven't already resolved this constraint
      # Find all items equal and complementary to x to resolve this constraint
      # Sets of all items eq and wc to x, we will build it iteratively
      eq_set = set(eq[x]);  eq_set.add(x)
      wc_set = set(wc[x])
      # Sets of items whose constraints we have already considered
      eq_done = set()
      wc_done = set()
      # While everything isn't done, propagate constraints
      while (eq_set - eq_done) or (wc_set - wc_done):
        for y in (eq_set - eq_done):
          assert y in keys and y not in eq_all, (x, y)
          eq_set.update(eq[y])  # If z == y == x, then z == x
          wc_set.update(wc[y])  # If z ~= y == x, then z ~= x
          eq_done.add(y)
        for y in (wc_set - wc_done):
          assert y in keys and y not in eq_all, (x, y)
          wc_set.update(eq[y])  # If z == y ~= x, then z ~= x
          eq_set.update(wc[y])  # If z ~= y ~= x, then z == x
          wc_done.add(y)
      # Store these complete constraints with each item
      for y in eq_set:
        eq_all[y] = eq_set
        wc_all[y] = wc_set
      for y in wc_set:
        eq_all[y] = wc_set
        wc_all[y] = eq_set
  
  return eq_all, wc_all

