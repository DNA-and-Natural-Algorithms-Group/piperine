import os
import subprocess
import unittest
import sys
from tempfile import mkstemp, mkdtemp
from .test_data import fixed_file
from time import time
import csv

from peppercompiler import compiler
from peppercompiler.design import PIL_parser
import stickydesign as sd

from .. import designer, tdm
from ..Srinivas2017 import translation, energetics

energyfuncs = energetics.energyfuncs()

# From a stackoverflow, 16571150
if sys.version_info >= (3,0):
    from io import StringIO
else:
    from StringIO import StringIO

class Capturing(list):
    def __enter__(self):
        self._stdout = sys.stdout
        sys.stdout = self._stringio = StringIO()
        return self
    def __exit__(self, *args):
        self.extend(self._stringio.getvalue().splitlines())
        sys.stdout = self._stdout

class TestCommandline(unittest.TestCase):
    def setUp(self):
        # CRN used in this test
        self.crn_rxn = 'A + B -> A + D\n'
        # Establish random temporary filenames
        self.tdir = mkdtemp(prefix='piperine_test')
        fid, self.basename = mkstemp(dir=self.tdir)
        os.close(fid)
        endings = ['.crn', '.fixed', '.pil', '.mfe', '_strands.txt', '{}.seqs']
        self.filenames = [ self.basename + suf for suf in endings ]
        self.crn, self.fixed, self.pil, self.mfe, self.strands, self.seqs = self.filenames
        fid, self.fixedscore = mkstemp(suffix='_fixed_score.csv', dir=self.tdir)
        os.close(fid)
        fid, self.reportfile = mkstemp(dir=self.tdir)
        os.close(fid)
        # Write CRN to basename.crn
        with open(self.crn, 'w') as f:
            f.write(self.crn_rxn)
        # Modules and module strings for import tests
        proc = subprocess.Popen(['piperine-design {} -n 3 -D -q'.format(self.crn)], stdout=subprocess.PIPE, shell=True)
        (out, err) = proc.communicate()
        self.ef = energetics.energyfuncs(targetdG=7.7)
        self.trans = translation

    def tearDown(self):
        filenames = os.listdir(self.tdir)
        for f in filenames:
            os.remove(os.path.join(self.tdir,f))
        os.rmdir(self.tdir)

    def runTest(self):
        self.test_run_designer_accepts_string_modules()

    def test_utility_access(self):
        command = ['piperine-design -h', 'piperine-score -h', 'piperine-select -h']
        for cmd in command:
            proc = subprocess.Popen([cmd], stdout=subprocess.PIPE, shell=True)
            (out, err) = proc.communicate()
            self.assertTrue('usage' in str(out), "Command failed: {}".format(cmd))

    def test_design_options(self):
        cmd_base = "piperine-design {} -D -q ".format(self.crn)
        cmd_options = ['-l 8 ', '-d 1', '-m 0.5', '-p 8 16 3']

        cmd_option = cmd_options[0]
        proc = subprocess.Popen([cmd_base + cmd_option], stdout=subprocess.PIPE, shell=True)
        (out, err) = proc.communicate()
        fixed_content = compiler.load_fixed(self.fixed)
        msg = "Commandline option -l {} makes toehold length {}. Fixed file:{}"
        self.assertEqual(8, len(fixed_content[0][-1]), msg.format(5, len(fixed_content[0][-1]), self.fixed))

        cmd_option = cmd_options[1]
        proc = subprocess.Popen([cmd_base + cmd_option], stdout=subprocess.PIPE, shell=True)
        (out, err) = proc.communicate()
        fixed_content = compiler.load_fixed(self.fixed)
        ths_strs = [''.join(tmp[-1]).lower() for tmp in fixed_content]
        th_sets  = [ths_strs[i:i+2] for i in range(0,len(ths_strs),2)]
        e_err, e_dev = energyfuncs.score_toeholds(th_sets)
        msg = "Commandline option {} makes energy deviation {}. Fixed file:{}".format(cmd_option, e_err, self.fixed)
        self.assertTrue(e_err < 1.0, msg)

        cmd_option = cmd_options[3]
        proc = subprocess.Popen([cmd_base + cmd_option], stdout=subprocess.PIPE, shell=True)
        (out, err) = proc.communicate()
        pil_content = PIL_parser.load_spec(self.pil)
        lengths = [8, 16-3, 3]
        domains = ['r0-toe-sa', 'r0-am', 'r0-cam']
        for length, domain in zip(lengths, domains):
            msg = "Commandline option {} makes domain {} with {} bases. PIL file:{}"
            msg = msg.format(cmd_option, domain, pil_content.base_seqs[domain].length, self.pil)
            self.assertEqual(length, pil_content.base_seqs[domain].length, msg)

    def test_score_options(self):
        cmd_base = "piperine-score {} ".format(self.crn)
        temp_str = ""
        for i in range(3):
            temp_str = temp_str + self.seqs.format(i) + " "
            proc = subprocess.Popen([cmd_base + temp_str + " -D -q -f {}".format(self.fixedscore)], stdout=subprocess.PIPE, shell=True)
            (out, err) = proc.communicate()
            self.assertTrue(os.path.exists(self.fixedscore))
            with open(self.fixedscore) as csvfile:
                reader = csv.reader(csvfile)
                nrows = 0
                for row in reader:
                    nrows += 1
                self.assertTrue(nrows, i+1)

    def test_select_options(self):
        fnames = os.listdir(self.tdir)
        csvs = [os.path.join(self.tdir, x) for x in fnames if '.csv' in x]
        cmd_base = "piperine-select "
        for csv_file in csvs:
            self.assertTrue(os.path.exists(csv_file))
            cmd_base = cmd_base + csv_file + " "
        cmd_base = cmd_base + " -f {} ".format(self.reportfile)
        proc = subprocess.Popen([cmd_base], stdout=subprocess.PIPE, shell=True)
        (out, err) = proc.communicate()
        self.assertTrue(os.path.exists(self.fixedscore))
        self.assertTrue(os.path.exists(self.reportfile))

def suite():
    tests = ['test_utility_access', 'test_design_options', 'test_score_options', 'test_select_options']
    # tests = ['test_select_options']
    return unittest.TestSuite(list(map(TestCommandline, tests)))
