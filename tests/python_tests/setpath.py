import sys

sys.path.append('/home/jmp/Dropbox/CRNtoDNASoftware/')
import crn_compiler.tests.run_unittests as ru
from crn_compiler.tests import TDMTests as tdmt

t = tdmt.TestTDM()
t.setUp()
t.runTest()
