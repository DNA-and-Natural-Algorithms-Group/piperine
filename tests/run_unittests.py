import unittest

import import_test
import TDMTests
import CompilationTests
import RunDesignerTest
import DSDClassesTests

def runem():
    suite = import_test.suite()
    unittest.TextTestRunner(verbosity=2).run(suite)

def run_tdm():
    suite = TDMTests.suite()
    unittest.TextTestRunner(verbosity=2).run(suite)

def run_compilation():
    suite = CompilationTests.suite()
    unittest.TextTestRunner(verbosity=2).run(suite)

def run_designer():
    suite = RunDesignerTest.suite()
    unittest.TextTestRunner(verbosity=2).run(suite)

def run_translation():
    suite = DSDClassesTests.suite()
    unittest.TextTestRunner(verbosity=2).run(suite)

def run_all():
    alltests = unittest.TestSuite(
        [
            x.suite() for x in 
                [import_test, TDMTests, CompilationTests, RunDesignerTest]
        ])
    unittest.TextTestRunner(verbosity=2).run(alltests)
