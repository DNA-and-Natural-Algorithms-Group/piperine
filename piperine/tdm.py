from __future__ import division, print_function

import os
import sys
from tempfile import mkstemp, mkdtemp

nupackpath = os.environ['NUPACKHOME']+'/bin/'

import numpy as np
import re
import random

from .designer import default_energyfuncs
whiteSpaceSearch = re.compile('\s+')

class MyProgress(object):
    class ImproperInput(Exception):
        pass

    def __init__(self, clockmax):
        self.clock = 0
        self.count = 0
        self.clockmax = clockmax
        self.update = 0.1 * self.clockmax
        if clockmax < 1:
            print('Input to constructor < 1')
            raise self.ImproperInput()
        print('Begin')

    def inc(self):
        self.clock += 1
        self.count += 1
        if self.count > self.update:
            prcnt = self.clock / self.clockmax
            print('Progress : {:.2%}'.format(prcnt))
            self.count = 0
        if self.clock >= self.clockmax:
            print('DONE')

def read_design(filename):
  from peppercompiler.nupack_out_grammar import document
  """Extracts the designed sequences and the mfe structures"""
  if not os.path.isfile(filename):
    error("Cannot load design. No such file '%s'." % filename)
  stats, total_n_star = document.parseFile(filename)
  seqs = {}
  structs = {}
  # Catches everything (struct, seq and seq*) except bad inters (struct N struct)
  #print stats
  for stat in stats:
    #print stat
    name, seq, n_star, gc, mfe, ideal_struct, actual_struct = stat
    seqs[name] = seq
    structs[name] = ideal_struct
  return [seqs, structs]

def Read_Finished(filename):
    f = open(filename,'r')
    lines = f.readlines()
    f.close()

    sequences = {}
    strands = {}
    for line in lines:
        # Check for comment
        if '#' in line:
            line = line.split('#')[0]
        parts = line.split(' ')
        if len(parts) == 4:
            name = parts[1]
            namseq = parts[3]
            namlen = len(namseq)
            if parts[0] == 'sequence':
                sequences[name] = namseq[0:namlen-1]
            elif parts[0] == 'strand':
                strands[name] = namseq[0:namlen-1]
    return (sequences, strands)

def make_pepper_seq_dict(pepperlist, seq_dict, update=False):
    Strandseq_dict = {}
    s_keys = list(seq_dict.keys())
    s_keys.sort(reverse=True)
    for pname in pepperlist:
        pswap = pname[:]
        pmem = pname[:]
        for key in s_keys:
            if key in pname and 'N' not in seq_dict[key]:
                pswap = pswap.replace(key, seq_dict[key])
                pmem = pmem.replace(key, '')
        if not update or pmem == '':
            Strandseq_dict[pname] = pswap

    return Strandseq_dict

def compare_sequence_notoe(s1, s2, toeholds):
     siz = min((len(s1), len(s2)));
     maxmatchsize = 0;
     ismaxmatch = "FALSE";
     mm_i = -1
     mm_j = -1
     for ll in range(5, siz + 1):
         for i in range(0, len(s1) - ll + 1):
             substr1 = s1[i:i+ll];
             for j in range(0, len(s2) - ll + 1):
                 substr2 = s2[j:j+ll];
                 if ((substr1==substr2) and substr1 not in toeholds):
                     maxmatchsize = ll;
                     mm_i = i
                     mm_j = j
                     if (maxmatchsize==siz):
                         ismaxmatch = "TRUE";
     return [ismaxmatch, maxmatchsize, mm_i, mm_j]

def get_seq_dicts(basename, heuristics_inputs, mfe_file=None, seq_file=None):
    if mfe_file is None:
        mfe_file = basename + ".mfe"
    if seq_file is None:
        seq_file = basename + ".seq"
    pep_sequences, pep_strands = Read_Finished(seq_file)
    domains_list = list(pep_sequences.keys())
    mfe_dict, cmplx_dict = read_design(mfe_file)
    seq_dict = mfe_dict.copy()
    seq_dict.update(pep_sequences)
    seq_dict.update(pep_strands)

    TopStrandlist, complex_names, BaseStrandlist, TopStranddict, BMlist, \
    NotToInteract = heuristics_inputs

    allpepperlist = list()
    for l in [TopStrandlist, complex_names, list(seq_dict.keys()), BMlist]:
        allpepperlist.extend(l)

    for l in list(NotToInteract.values()):
        allpepperlist.extend(l)

    allpepperlist.extend([ s[:-1] for s in BaseStrandlist] )

    seq_dict = make_pepper_seq_dict(allpepperlist, seq_dict, update=True)

    return seq_dict, cmplx_dict, domains_list

def get_heuristics_inputs(gates, strands):
    # Generate pepper-lists
    TopStrandlist = []
    complex_names = []
    TopStranddict = dict()
    BMlist = []
    NotToInteract = dict()
    for strand in strands:
        TopStrandlist.extend(strand.get_top_strands())
        BMlist.extend(strand.get_bms())

    BaseStrandlist = []
    for gate in gates:
        bases = gate.get_base_domains()
        for th in bases:
            th_c = th+'*'
            if th_c not in BaseStrandlist:
                BaseStrandlist.append(th_c)
        TopStrandlist.extend(gate.get_top_strands())
        complex_names.extend(gate.get_complexes())
        TopStranddict.update(gate.get_top_strand_dict())

    for gate in gates:
        for th_c in BaseStrandlist:
            if th_c not in NotToInteract:
                NotToInteract[th_c] = gate.get_noninteracting_peppernames(th_c[:-1])
            else:
                NotToInteract[th_c].extend(gate.get_noninteracting_peppernames(th_c[:-1]))

    for strand in strands:
        for th_c in BaseStrandlist:
            if th_c not in NotToInteract:
                NotToInteract[th_c] = strand.get_noninteracting_peppernames(th_c[:-1])
            else:
                NotToInteract[th_c].extend(strand.get_noninteracting_peppernames(th_c[:-1]))

    return (TopStrandlist, complex_names, BaseStrandlist, TopStranddict, BMlist,
            NotToInteract)

def EvalCurrent(basename, gates, strands, compile_params=(7, 15, 2),
                header=True, testname=None, seq_file=None, mfe_file=None,
                quick=False, energyfuncs=default_energyfuncs,
                includes=None, clean=True):
    if not testname:
        testname = basename
    if not seq_file:
        seq_file = basename + '.seq'
    if not mfe_file:
        mfe_file = basename + '.mfe'

    heuristics_inputs = get_heuristics_inputs(gates, strands)
    (TopStrandlist, complex_names, BaseStrandlist, TopStranddict, BMlist,
     NotToInteract) = heuristics_inputs

    seq_dict, cmplx_dict, domains_list = get_seq_dicts(basename, heuristics_inputs,
                                                       mfe_file, seq_file)


    print('Start WSI computation')
    if quick:
        ssm_scores = np.random.rand(6)
    else:
        ssm_scores = Spurious_Weighted_Score(basename, domains_list, seq_dict,
                                             compile_params=compile_params,
                                             includes=includes, clean=clean)
    ssm_names = ['WSAS', 'WSIS', \
                 'WSAS-M', 'WSIS-M', \
                 'Verboten', 'Spurious']
    print('')

    # Score cross-strand spurious interactions
    print('Start Cross-Strand spurious interactions computation')
    if quick:
        css_scores = np.random.rand(4)
    else:
        css_scores  = NUPACK_Eval(seq_dict, TopStrandlist, BaseStrandlist, \
            NotToInteract, ComplexSize = 2, T = 25.0, material = 'dna',\
             clean=clean, quiet=True)
    css_names = ['TSI avg', 'TSI max', \
                 'TO avg', 'TO max']
    print('')

    print('Start bad nucleotide percent computation')
    if quick:
        ted_scores = [np.random.rand(), 'BAD', np.random.rand()]
    else:
        ted_scores = NUPACK_Eval_bad_nucleotide(seq_dict, cmplx_dict, complex_names,\
            prefix='tube_ensemble', clean=clean)
    ted_names = ['BN% max', 'Max Defect Component',
                 'BN% avg']
    print('')

    # Retrieve toeholds for BM score calculation
    th_strs = [ s.get_ths() for s in strands ]
    toeholds = [seq_dict[i] for ths in th_strs for i in ths]

    # Score weighted spurious interactions
    print('Start BM score computation')
    if quick:
        bm_scores = np.random.rand(2)
    else:
        bm_scores = BM_Eval(seq_dict, BMlist, toeholds)
    bm_names = ['WS-BM', 'Max-BM']
    print('')

    # Score intra-strand spurious interactions and toehold availability
    print('Start Single-Strand spurious score computation')
    if quick:
        ss_scores = np.random.rand(4)
    else:
        ss_scores = SS_Eval(seq_dict, TopStranddict, T = 25.0, material = 'dna', clean=clean)
    ss_names = ['SSU min', 'SSU avg', 'SSTU min', 'SSTU avg']
    print('')

    # Convert the list of toeholds from each strand into nucleotide sequences
    toeholds = [ [seq_dict[i] for i in s.get_ths()] for s in strands ]

    if quick:
        th_scores = np.random.random((2,))
    else:
        th_scores = energyfuncs.score_toeholds(toeholds)
    th_names = ['dG Error', 'dG Range']

    score_list = [css_scores, bm_scores, ss_scores, ted_scores,
                  ssm_scores, th_scores]
    names_list = [css_names, bm_names, ss_names, ted_names,
                  ssm_names, th_names]
    scores = [ elem for sub in score_list for elem in sub]
    names  = [ elem for sub in names_list for elem in sub]

    if header:
        output = (scores, names)
    else:
        output = scores
    return output

def NUPACK_Eval_bad_nucleotide(mfe_seqs, ideal_structs, complex_names, \
                            ComplexSize=3, T=25.0, material='dna', \
                            clean=True, quiet=True, prefix='ted_calc',
                            tmpdir=None):
    # This function takes in the design's sequences and intended interaction
    # structures and returns the tube ensemble defect, the concentration of incorr
    # ectly base-paired nucleotides in a tube. I think this is best prepared and t
    # hought of not as a single tube, but as if each intermediate were prepared in
    #  separate tubes and the concentration is summed over each of these tubes.

    # Loop through the complexes in all systems to calculate their contribution
    # to the tube ensemble defect (ted) , the running sum of such contributions
    # Set up empty lists
    ted_vec = np.empty(len(complex_names))
    name_list = list()
    bp_vec = np.empty(len(complex_names))
    # Set parameters
    target_conc = 1e-06
    counter = 0

    prog = MyProgress(len(complex_names))
    for cmpx in complex_names:
        # Grab complex sequences and structure
        cmpx_name = cmpx
        seq_list = mfe_seqs[cmpx_name].split('+')
        seq = mfe_seqs[cmpx_name].replace('+', '')
        nseq = len(seq_list)
        params = [nseq, T, material, quiet, prefix]
        struct = ideal_structs[cmpx_name]
        bp = len(re.findall('[()]', struct))
        # Retrieve the estimated complex concentration
        est_conc = NUPACK_Cmpx_Conc(seq_list, params, clean=clean, tmpdir=tmpdir)
        cmpx_defect = NUPACK_Cmpx_Defect(seq_list, struct, params, clean=clean, tmpdir=tmpdir)
        ted = cmpx_defect * min(est_conc, target_conc) + \
               len(seq) * max(target_conc - est_conc, 0)
        ted_vec[counter] = ted
        name_list.append(cmpx_name)
        bp_vec[counter] = bp
        counter += 1
        prog.inc()

    # Convert TED to % bad nucleotides out of total nucleotides
    bad_nuc_vec = 100 * ted_vec / (target_conc * bp_vec)
    bad_nuc_max = bad_nuc_vec.max()
    bn_max_name = name_list[int(np.where(bad_nuc_vec == bad_nuc_max)[0][0])]
    #return [Bad Nucleotide % max, max complex name, mean bad nuc]
    return [bad_nuc_max, bn_max_name, bad_nuc_vec.mean()]

def NUPACK_Eval(seq_dict, TopStrandlist, BaseStrandlist, NotToInteract,\
                ComplexSize = 2, T = 25.0, material = 'dna', \
                clean=True, quiet=True):
    numstrands = len(TopStrandlist)

    TopSpuriousPairwise = np.zeros([numstrands, numstrands]);

    print('Calculating Top Strand Pairwise interactions')
    countmax = (numstrands**2 + numstrands)/2
    prog = MyProgress(countmax)
    for i in range(numstrands):
        for j in range(i,numstrands):
            intij = NUPACKIntScore(TopStrandlist[i], TopStrandlist[j],
                                   seq_dict, ComplexSize, T, material, quiet, clean=clean)
            TopSpuriousPairwise[i, j] = intij
            TopSpuriousPairwise[j, i] = intij
            prog.inc()

    numbase = len(BaseStrandlist)

    BaseSpurious = np.zeros([numbase, 1]);

    print('Calculating Toehold occupation')
    countmax = sum(map(len, list(NotToInteract.values())))
    prog = MyProgress(countmax)
    for i in range(numbase):
        thisstrand = BaseStrandlist[i]
        notinteract = NotToInteract[thisstrand]
        nonintnum = len(notinteract)
        for j in range(nonintnum):
            intij = NUPACKIntScore(thisstrand, notinteract[j], seq_dict,
                                   ComplexSize, T, material, quiet, clean=clean)
            BaseSpurious[i] = BaseSpurious[i] + intij
            prog.inc()

    TSI_vec = TopSpuriousPairwise.sum(0)
    return [TSI_vec.mean(), TSI_vec.max(), BaseSpurious.mean(),
            BaseSpurious.max()]


def SS_Eval(seq_dict, TopStranddict, T = 25.0, material = 'dna', clean=True):
    numstrands = len(TopStranddict)

    MinProbs = []
    avg_Unpaired = []

    MinProbs_toe = []
    avg_Unpaired_toe = []

    prog = MyProgress(numstrands)
    for strand, toe_regions in TopStranddict.items():
            [min_Unpaired, sum_Unpaired, min_Unpaired_toe, sum_Unpaired_toe,\
             NumBases] = \
                NUPACKSSScore(strand, seq_dict, T, material, toe_regions, clean=clean)
            # Grow the minimum unpaired probability lists
            MinProbs.append(min_Unpaired)
            MinProbs_toe.append(min_Unpaired_toe)

            # Grow the average strand unpaired probability lists
            avg_Unpaired.append(sum_Unpaired/NumBases)
            if len(toe_regions) > 0:
                avg_Unpaired_toe.append(sum_Unpaired_toe/len(toe_regions))
            else:
                avg_Unpaired_toe.append(0)
            prog.inc()


    return [min(MinProbs), np.mean(avg_Unpaired), \
            min(MinProbs_toe), np.mean(avg_Unpaired_toe)]

def BM_Eval(seq_dict, BMlist, toeholds):
    w_exp = np.concatenate([np.zeros((5,)), np.power(2, np.arange(6))])
    BM_score = 0
    Largest_match = 0

    numstrings = len(BMlist);

    prog = MyProgress((numstrings**2 - numstrings)/2)
    for ctr in range(numstrings):
        strand1 = BMlist[ctr];
        for strand2 in BMlist[ctr+1:]:
            [ismaxmatch, maxmatch, mm_i, mm_j] = \
                compare_sequence_notoe(seq_dict[strand1],
                                       seq_dict[strand2],
                                       toeholds)
            if maxmatch > Largest_match:
                Largest_match = maxmatch

            BM_score = BM_score + w_exp[int(min(maxmatch, 10))]
            prog.inc()

    return [BM_score, Largest_match]

def Spurious_Weighted_Score(basename,
                            domains_list,
                            seq_dict,
                            compile_params=(7, 15, 2),
                            bored=6,
                            tmax=15,
                            spurious_range=10,
                            beta=5,
                            clean=True,
                            includes=None,
                            tmpdir=None):
    # This function calculates Niranjan's weighted spurious interaction score, for i
    # nteractions between k and km nucleotides long. It takes the following steps:
    #   * Read in sequences from mfe file
    #   * Create fixed file for all sequences except clamps
    #   * Recompile circuit
    #   * Run spuriousSSM negative design and generate scores
    #   * Grab spuriousSSM scores and calculate NSIH

    import sys
    import subprocess
    import re
    import os
    import numpy as np

    from .designer import call_compiler, call_design

    # Command parameters
    if tmpdir is None:
        fid, prefix = mkstemp()
    else:
        tdir = mkdtemp(dir=tmpdir)
        fid, prefix = mkstemp(dir=tdir)
    ssm_params = "bored=%s tmax=%s spurious_range=%s" % (bored, tmax, spurious_range)

    fid, fixed_file = mkstemp(suffix='.fixed', dir=tmpdir)
    os.close(fid)
    fid, compiled_file = mkstemp(suffix='.pil', dir=tmpdir)
    os.close(fid)
    fid, save_file = mkstemp(suffix='.save', dir=tmpdir)
    os.close(fid)
    fid, out_file = mkstemp(suffix='.mfe', dir=tmpdir)
    os.close(fid)
    fid, spurious_output = mkstemp(suffix='.txt', dir=tmpdir)
    os.close(fid)

    # Write sequences to fixed file
    f = open(fixed_file, 'w')
    for dom in domains_list:
        f.write('sequence ' + dom + ' = ' + seq_dict[dom] + ' \n')
    f.close()

    # Compile to a pil file
    call_compiler(basename,
                    args=compile_params,
                    outputname=compiled_file,
                    savename=save_file,
                    fixed_file=fixed_file,
                    includes=includes)

    # Make constraint files
    fid, design_tmp = mkstemp(dir=tmpdir)
    os.close(fid)
    call_design(basename, infilename=compiled_file, outfilename=out_file,
                just_files=True, tempname=design_tmp)
    stname = design_tmp+'.st'
    wcname = design_tmp+'.wc'
    eqname = design_tmp+'.eq'

    # Commented code generates MFE file
    spur_exc = 'spuriousSSM score=automatic template=%s wc=%s eq=%s %s> %s'
    command = spur_exc % (stname, wcname, eqname, ssm_params, spurious_output)
    os.system(command)
    #subprocess.check_call(spur_exc)

    w_lin = np.concatenate([np.zeros((beta, )), np.arange(12-beta)+1, 7*np.ones((2,))])

    f = open(stname)
    st = f.readline()[:-1]
    f.close()
    num_strands = len(re.findall('\s+', st))+1

    # re string for decimal/integer scraping
    num = "\d+\.?\d*"

    # Read in entire file, search for instances of "spurious1"
    spc_text = open(spurious_output).read()
    lines = spc_text[spc_text.rfind('spurious1'):].split('\n')
    i, j = [ int(x) for x in re.findall(num, lines[0])[1:]]
    w = w_lin[i-1:j]
    vec = np.array([ np.float(x) for x in re.findall(num, lines[2])])
    score_vec = vec * w
    mis_intra_score = score_vec.sum()  #/ num_strands

    vec = np.array([ np.float(x) for x in re.findall(num, lines[4])])
    score_vec = vec * w
    mis_inter_score = score_vec.sum()  #/ num_strands

    lines = spc_text[spc_text.rfind('spurious('):].split('\n')
    i, j = [ int(x) for x in re.findall(num, lines[0])]
    w = w_lin[i-1:j]
    vec = np.array([ np.float(x) for x in re.findall(num, lines[2])])
    score_vec = vec * w
    spc_intra_score = score_vec.sum()  #/ num_strands

    vec = np.array([ np.float(x) for x in re.findall(num, lines[4])])
    score_vec = vec * w
    spc_inter_score = score_vec.sum()  #/ num_strands

    vec_str = spc_text[spc_text.rfind('** score_verboten'):]
    verboten_score = np.float(re.findall(num, vec_str)[0])  #/ num_strands

    vec_str = spc_text[spc_text.rfind('-weighted score = '):].split('\n')[0]
    wsi_score = np.float(re.findall(num, vec_str)[-1])  #/ num_strands

    if clean:
        for f in [fixed_file, compiled_file, save_file, out_file,
                  stname, wcname, eqname, spurious_output]:
            os.remove(f)

    return [spc_intra_score, spc_inter_score, \
            mis_intra_score, mis_inter_score, \
            verboten_score, wsi_score]

def NUPACK_Cmpx_Conc(seqs, params=[3, 25, 'dna', 1, 'ted_calc'], clean=True, tmpdir=None):
    # This function calls the NUPACK methods 'complexes' and 'concentrations' to
    # calculate the expected concentration of the desired complex. Note that the
    # output of 'concentrations' will contain the expected concentrations of all
    # possible complexes, but this function returns only that of the desired com
    # plex, which is always 1 1 1.
    # Retrieve parameters
    ComplexSize, T, material, quiet, f_prefix = params
    # Setup input files to nupack commands
    intsc = 0;
    if tmpdir is None:
        fid, prefix = mkstemp()
    else:
        tdir = mkdtemp(dir=tmpdir)
        fid, prefix = mkstemp(dir=tdir)
    os.close(fid)
    ofile = prefix + '.out'
    cfile = prefix + '.con'
    ifile = prefix + '.in'
    efile = prefix + '.eq'
    xfile = prefix + '.osx'
    file_list = [ofile, cfile, ifile, efile, xfile, prefix]
    #for fn in file_list:
    #    with open(fn, 'w') as f:
    #        f.write('')
    #        f.close()

    # Complexes input file
    n_seqs = len(seqs)
    f = open(ifile,'w')
    f.write(str(n_seqs) + '\n')
    for seq in seqs:
      f.write(seq + '\n')

    f.write(str(ComplexSize))
    f.close()

    # Concentrations input file, start with 1uM for all
    f = open(cfile,'w')
    for i in range(n_seqs):
      f.write('1e-6\n')

    f.close()

    # Run NUPACK's complexes function
    cmd = nupackpath + 'complexes -T %.1f '+\
          '-material %s -ordered -pairs -mfe -dangles some -sodium 0.5 %s > %s'
    cmd = cmd % (T,material,prefix,ofile)
    os.system(cmd)

    # Run NUPACK's concentrations function
    if quiet:
        cmd = nupackpath + 'concentrations '+\
              '-ordered -pairs -cutoff 0.000001 -sort 0 -quiet ' + prefix
    else:
        cmd = nupackpath + 'concentrations '+\
              '-ordered -pairs -cutoff 0.000001 -sort 0 ' + prefix
    os.system(cmd)

    # Read the equilibrium concentrations
    f = open(efile)
    lines = f.readlines()
    f.close()

    if clean:
        for f in file_list:
            if os.path.isfile(f):
                os.remove(f)

    # Count the number of commented lines
    ctr = 0;
    while (lines[ctr][0] == '%'):
        ctr = ctr + 1

    # Scrape the concentrations from the output file
    out = 0.0
    for count in range(ctr, len(lines)):
        lineData = whiteSpaceSearch.split(lines[count])
        if lineData[2:(ComplexSize+2)] == (ComplexSize ) * ['1'] :
            trial_out = float(lineData[-2])
            if trial_out > out:
                out = trial_out

    return out


def NUPACK_Cmpx_Defect(seqs, struct, params=[3, 25, 'dna', 1, 'ted_calc'], clean=True,
                       tmpdir=None):
    '''
    This function calls the NUPACK methods 'defects' and reports the average num
    ber of basepairs that do not match their intended state.
    '''
    ## TODO : Add error for when complex size is less than the number of provided sequences
    # Retrieve parameters
    ComplexSize, T, material, quiet, prefix = params
    cmpx_string = ' '.join(["{}".format(x+1) for x in range(ComplexSize)]) + '\n'
    # Setup input files to nupack commands
    intsc = 0;
    if tmpdir is None:
        fid, prefix = mkstemp()
    else:
        tdir = mkdtemp(dir=tmpdir)
        fid, prefix = mkstemp(dir=tdir)
    os.close(fid)
    ofile = prefix + '.out'
    ifile = prefix + '.in'
    file_list = [prefix, ofile, ifile]

    # Complexes input file
    n_seqs = len(seqs)
    f = open(ifile,'w')
    f.write(str(n_seqs) + '\n')
    for seq in seqs:
      f.write(seq + '\n')

    # f.write('1 2 3\n')
    f.write(cmpx_string)
    f.write(struct)
    f.close()

    # Run NUPACK's defect function
    cmd = nupackpath + 'defect -T %.1f '+\
          '-material %s -dangles some -multi -sodium 0.5 %s > %s'
    cmd = cmd % (T,material,prefix,ofile)
    os.system(cmd)

    # Read the equilibrium concentrations
    f = open(ofile)
    lines = f.readlines()
    f.close()

    if clean:
        for f in file_list:
            if os.path.isfile(f):
                os.remove(f)

    # Count the number of commented lines
    ctr = 0
    while (lines[ctr][0] == '%'):
        ctr = ctr + 1
    return float(lines[ctr])


def NUPACKIntScore(str1, str2, seq_dict,
                   ComplexSize=2,
                   T=25.0,
                   material='dna',
                   quiet=True,
                   clean=True,
                   tmpdir=None):
    # Formerly IntScore
    import os

    intsc = 0;
    if tmpdir is None:
        fid, fname = mkstemp()
    else:
        tdir = mkdtemp(dir=tmpdir)
        fid, fname = mkstemp(dir=tdir)
    os.close(fid)
    ofile = fname + '.out'
    cfile = fname + '.con'
    ifile = fname + '.in'
    efile = fname + '.eq'
    file_list = [fname, ofile, cfile, ifile, efile]

    f = open(ifile,'w')
    f.write('2\n')
    f.write(seq_dict[str1] + '\n')
    f.write(seq_dict[str2] + '\n')
    f.write(str(ComplexSize))
    f.close()

    f = open(cfile,'w')
    f.write('1e-6\n')
    f.write('1e-6\n')
    f.close()

    cmd = nupackpath + 'complexes -T %.1f '+\
          '-material %s -ordered -pairs -mfe -dangles some -sodium 0.5 %s > %s'
    cmd = cmd % (T,material,fname,ofile)
    os.system(cmd)
    #cmd = '/home/juic/Downloads/PROGS/nupack3.0.4/bin/concentrations
    # -ordered -pairs -cutoff 0.000001 -sort 0 -quiet ' + fname
    if quiet:
        cmd = nupackpath + 'concentrations '+\
              '-ordered -pairs -cutoff 0.000001 -sort 0 -quiet ' + fname
    else:
        cmd = nupackpath + 'concentrations '+\
              '-ordered -pairs -cutoff 0.000001 -sort 0 ' + fname
    os.system(cmd)

    f = open(efile)
    lines = f.readlines()
    f.close()
    ctr = 0;
    while (lines[ctr][0] == '%'):
        ctr = ctr + 1
    ctr = ctr + 2;

    if clean:
        for f in file_list:
            if os.path.isfile(f):
                os.remove(f)

    for count in range(ctr, len(lines)):
        lineData = whiteSpaceSearch.split(lines[count])
        rawsc = float(lineData[-2])
        intsc = intsc + 100*rawsc/(2*1e-6)

    return intsc

def NUPACKSSScore(str1, seq_dict, T=25.0, material='dna', toe_region=[None], clean=True, tmpdir=None):
    # This function calls NUPACK's complexes and concentrations functions, which
    # allows one to calculate the unpaired probabilities of each nucleotide across
    # the most common secondary structures.
    import os
    ComplexSize = 1
    min_Unpaired = 1;
    sum_Unpaired = 0;

    min_Unpaired_toe = 1;
    sum_Unpaired_toe = 0;

    # Grab the sequence to assay
    seq = seq_dict[str1]
    UnpairedIndex = len(seq) + 1

    # Prepare filenames
    if tmpdir is None:
        fid, fname = mkstemp()
    else:
        tdir = mkdtemp(dir=tmpdir)
        fid, fname = mkstemp(dir=tdir)
    os.close(fid)
    ofile = fname + '.out'
    ofile2 = fname + '.con.out'
    cfile = fname + '.con'
    ifile = fname + '.in'
    pfile = fname + '.fpairs'
    file_list = (fname, ofile, ofile2, cfile, ifile, pfile)

    # Write the number of strands, the sequence, and the maximum complex size
    # to the input file, to be read by NUPACK
    f = open(ifile,'w')
    f.write('1\n')
    f.write(seq + '\n')
    f.write(str(ComplexSize))
    f.close()

    f = open(cfile,'w')
    f.write('1e-6\n')
    f.close()

    cmd = nupackpath + 'complexes -T %.1f '+\
          '-material %s -ordered -pairs -mfe -dangles some -sodium 0.5 %s > %s'
    cmd = cmd % (T,material,fname,ofile)
    os.system(cmd)
    cmd = nupackpath + 'concentrations  '+\
          '-ordered -pairs -cutoff 0.000001 -sort 0 -quiet ' + fname
    os.system(cmd)

    # Read in nucleotide pairing information
    f = open(pfile)
    lines = f.readlines()
    f.close()

    if clean:
        for f in file_list:
            if os.path.isfile(f):
                os.remove(f)

    ctr = 0
    while (lines[ctr][0] == '%'):
        ctr = ctr + 1
    ctr = ctr + 1

    for count in range(ctr, len(lines)):
        lineData = whiteSpaceSearch.split(lines[count])
        if lineData[1]==str(UnpairedIndex):
            frac_unp = float(lineData[-2])
            sum_Unpaired = sum_Unpaired + frac_unp
            if frac_unp < min_Unpaired:
                min_Unpaired = frac_unp
            if int(lineData[0]) in toe_region:
                sum_Unpaired_toe = sum_Unpaired_toe + frac_unp
                if frac_unp < min_Unpaired_toe:
                    min_Unpaired_toe = frac_unp

    #avg_Unpaired = float(sum_Unpaired)/(UnpairedIndex-1)
    if None in toe_region:
        return [min_Unpaired, sum_Unpaired, UnpairedIndex-1]
    else:
        return [min_Unpaired, sum_Unpaired, min_Unpaired_toe, \
                sum_Unpaired_toe, UnpairedIndex-1]
