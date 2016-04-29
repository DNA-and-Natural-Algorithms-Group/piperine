import import_test
import TDMTests
import CompilationTests
import unittest

def runem():
    suite = import_test.suite()
    unittest.TextTestRunner(verbosity=2).run(suite)

def run_tdm():
    suite = TDMTests.suite()
    unittest.TextTestRunner(verbosity=2).run(suite)

def run_compilation():
    suite = CompilationTests.suite()
    unittest.TextTestRunner(verbosity=2).run(suite)
