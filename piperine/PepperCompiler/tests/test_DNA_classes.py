#! /usr/bin/env python

import unittest

import DNA_classes

class TestDNAClasses(unittest.TestCase):
  
  ## Sequence tests
  example_sequence = ( 
    ( ("name", "prefix-", [(5, "N")], None), ("NNNNN", 5, False) ),
    ( ("name", "prefix-", [(5, "N")], 5), ("NNNNN", 5, False) ),
    ( ("name", "prefix-", [(1, "S"), (1, "W"), (1, "A"), (1, "C"), (1, "Y"), (1, "B")], None), 
        ("SWACYB", 6, False) ),
    ( ("name", "prefix-", [("?", "N")], 7), ("NNNNNNN", 7, False) ),
    ( ("name", "prefix-", [(2, "S"), ("?", "N"), (2, "S")], 7), ("SSNNNSS", 7, False) ),
    ( ("name", "prefix-", [], None), ("", 0, True) ),
  )
  
  def test_sequence(self):
    """Just check that we can make a sequence without error."""
    name, prefix, constraints, length = "name", "prefix-", [(5, "N")], None
    sequence = DNA_classes.Sequence(name, prefix, constraints, length)
    self.assert_(isinstance(sequence, DNA_classes.Sequence))
  
  def test_sequence_correct(self):
    """Check that we can make a sequence correctly."""
    name, prefix, constraints, length = "name", "prefix-", [(5, "N")], 5
    sequence = DNA_classes.Sequence(name, prefix, constraints, length)
    
    self.assert_(isinstance(sequence, DNA_classes.Sequence))
    self.assertEqual(name, sequence.name)
    self.assertEqual(prefix + name, sequence.full_name)
    self.assertEqual("NNNNN", sequence.const)
    self.assertEqual(length, sequence.length)
    self.assertEqual(False, sequence.reversed)
    self.assertEqual(False, sequence.dummy)
    
    self.assert_(isinstance(sequence.wc, DNA_classes.ReverseSequence))
    self.assertEqual(name + "*", sequence.wc.name)
    self.assertEqual(prefix + name + "*", sequence.wc.full_name)
    self.assertEqual(length, sequence.wc.length)
    self.assertEqual(True, sequence.wc.reversed)
    self.assertEqual(False, sequence.wc.dummy)
    
    self.assert_(sequence is sequence.wc.wc)
  
  def test_sequence_example(self):
    """Test example sequences and check that constraints were handled correctly."""
    for (name, prefix, constraints, length), (const_out, length_out, dummy) in self.example_sequence:
      sequence = DNA_classes.Sequence(name, prefix, constraints, length)
      self.assertEqual(const_out, sequence.const)
      self.assertEqual(length_out, sequence.length)
      self.assertEqual(dummy, sequence.dummy)
  
  
  ## Super sequence tests
  def test_super_sequence(self):
    """Just check that we can make a super-sequence without error."""
    name, prefix, constraints, length = "name", "prefix-", [[(5, "N")]], None
    sequence = DNA_classes.SuperSequence(name, prefix, constraints, length)
    self.assert_(isinstance(sequence, DNA_classes.SuperSequence))
  
  def test_super_sequence_correct(self):
    """Check that we can make a super-sequence correctly."""
    name, prefix, constraints, length = "name", "prefix-", [[(5, "N")]], 5
    sequence = DNA_classes.SuperSequence(name, prefix, constraints, length)
    
    self.assert_(isinstance(sequence, DNA_classes.SuperSequence))
    self.assertEqual(name, sequence.name)
    self.assertEqual(prefix + name, sequence.full_name)
    self.assertEqual(length, sequence.length)
    self.assertEqual(False, sequence.reversed)
    self.assertEqual(False, sequence.dummy)
    
    self.assertEqual(1, len(sequence.seqs))
    self.assert_(isinstance(sequence.seqs[0], DNA_classes.AnonymousSequence))
    self.assertEqual(5, sequence.seqs[0].length)
    
    self.assert_(isinstance(sequence.wc, DNA_classes.ReverseSuperSequence))
    self.assertEqual(name + "*", sequence.wc.name)
    self.assertEqual(prefix + name + "*", sequence.wc.full_name)
    self.assertEqual(length, sequence.wc.length)
    self.assertEqual(True, sequence.wc.reversed)
    self.assertEqual(False, sequence.wc.dummy)
    
    self.assert_(sequence is sequence.wc.wc)
  
  # TODO: Test with actual subsequences
  
  
  ## Strand tests
  def test_strand(self):
    """Just check that we can make a strand without error."""
    name, prefix, constraints, length = "name", "prefix-", [[(5, "N")]], None
    sequence = DNA_classes.Strand(name, prefix, constraints, length)
    self.assert_(isinstance(sequence, DNA_classes.Strand))
  
  # TODO: Test with actual subsequences
  
  
  ## TODO: Structure tests
  ## TODO: Kinetics tests

if __name__ == '__main__':
  unittest.main()

