import unittest

from . import CRNImportTests
from . import TDMTests
from . import CompilationTests
from . import RunDesignerTests
from . import Srinivas2017Tests
from . import TDM_NUPACK_tests
from . import CommandlineTests

def run_CRN_import_test():
    suite = CRNImportTests.suite()
    unittest.TextTestRunner(verbosity=2).run(suite)

def run_tdm():
    suite = TDMTests.suite()
    unittest.TextTestRunner(verbosity=2).run(suite)

def run_compilation():
    suite = CompilationTests.suite()
    unittest.TextTestRunner(verbosity=2).run(suite)

def run_designer():
    suite = RunDesignerTests.suite()
    unittest.TextTestRunner(verbosity=2).run(suite)

def run_translation():
    suite = DSDClassesTests.suite()
    unittest.TextTestRunner(verbosity=2).run(suite)

def run_nupack():
    suite = TDM_NUPACK_tests.suite()
    unittest.TextTestRunner(verbosity=2).run(suite)

def run_commandline():
    suite = CommandlineTests.suite()
    unittest.TextTestRunner(verbosity=2).run(suite)

def run_all():
    alltests = unittest.TestSuite(
        [
            x.suite() for x in
                #[import_test, TDMTests, CompilationTests, RunDesignerTest, DSDClassesTests, TDM_NUPACK_tests, CommandlineTests]
                [CommandlineTests, CRNImportTests, TDMTests, CompilationTests, RunDesignerTests, Srinivas2017Tests, TDM_NUPACK_tests]
        ])
    unittest.TextTestRunner(verbosity=2).run(alltests)
