from __future__ import division

import os
import sys

nupackpath = os.environ['NUPACKHOME']+'/bin/'

from numpy import array, hstack
import numpy as np
from scipy import zeros
import re
import random
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
            print 'Input to constructor < 1'
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

################################################################################
#
# Functions to read Template and Domain files
# \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
################################################################################

def read_design(filename):
  from PepperCompiler import nupack_out_grammar as ngram
  """Extracts the designed sequences and the mfe structures"""
  if not os.path.isfile(filename):
    error("Cannot load design. No such file '%s'." % filename)
  stats, total_n_star = ngram.document.parseFile(filename)
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
    s_keys = seq_dict.keys()
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

################################################################################
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# Functions to read Template and Domain files
# 
################################################################################

################################################################################
# 
# Sequence manipulation and mutation functions
# \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
################################################################################

def Complement(base):
    if base=='A':
        return 'T'
    elif base=='T':
        return 'A'
    elif base=='C':
        return 'G'
    elif base=='G':
        return 'C'

def Random_base(Type):
    if Type == 'W':
        return random.choice(['A', 'T'])
    elif Type == 'N':
        return random.choice(['A','C','T'])

def GetCurrentSub():
    current_Sub ={
    'W1':'G','W|1':'C',
    'W2':'C','W|2':'G',
    'W3':'G','W|3':'C',
    'N1':'C','N|1':'G',
    'N2':'G','N|2':'C',
    'N3':'G','N|3':'C',
    'N4':'G','N|4':'C',
    'N5':'G','N|5':'C',
    'N6':'G','N|6':'C',
    'N7':'G','N|7':'C',
    'N8':'A','N|8':'T',
    'N9':'G','N|9':'C',
     }
    return current_Sub

def GetRandomStart():
    W1 = Random_base(Type = 'W')
    W2 = Random_base(Type = 'W')
    N1 = Random_base(Type = 'N')
    N2 = Random_base(Type = 'N')
    N3 = Random_base(Type = 'N')
    N4 = Random_base(Type = 'N')
    N5 = Random_base(Type = 'N')
    N6 = Random_base(Type = 'N')
    return [W1, W2, N1, N2, N3, N4, N5, N6]

def GetBigRandomStart():
    W1 = Random_base(Type = 'W')
    W2 = Random_base(Type = 'W')
    W3 = Random_base(Type = 'W')
    N1 = Random_base(Type = 'N')
    N2 = Random_base(Type = 'N')
    N3 = Random_base(Type = 'N')
    N4 = Random_base(Type = 'N')
    N5 = Random_base(Type = 'N')
    N6 = Random_base(Type = 'N')
    N7 = Random_base(Type = 'N')
    return [W1, W2, W3, N1, N2, N3, N4, N5, N6, N7]


def MakeSub(W1, W2, W3, N1, N2, N3, N4, N5, N6, N7, N8, N9):
    tempSub = {
    'W1': W1,'W|1':Complement(W1),
    'W2': W2,'W|2':Complement(W2),
    'W3': W3,'W|3':Complement(W3),
    'N1': N1,'N|1':Complement(N1),
    'N2': N2,'N|2':Complement(N2),
    'N3': N3,'N|3':Complement(N3),
    'N4': N4,'N|4':Complement(N4),
    'N5': N5,'N|5':Complement(N5),
    'N6': N6,'N|6':Complement(N6),
    'N7': N7,'N|7':Complement(N7),
    'N8': N8,'N|8':Complement(N8),
    'N9': N9,'N|9':Complement(N9),
    }
    return tempSub

def MakeHPSub(X1, X2, X3, X4, X5):
    tempSub = {
    'X1': X1,'Xp1':Complement(X1),
    'X2': X2,'Xp2':Complement(X2),
    'X3': X3,'Xp3':Complement(X3),
    'X4': X4,'Xp4':Complement(X4),
    'X5': X5,'Xp5':Complement(X5)
    }
    return tempSub

def Mutate(InputBaseList):
    # Input List = [W1, W2, N1, N2, N3, N4, N5, N6]
    rand_index = random.randint(0, 7)
    BaseToMutate = InputBaseList[rand_index]

    if rand_index < 2:
        if BaseToMutate=='A':
            NewBase = 'T'
        elif BaseToMutate=='T':
            NewBase = 'A'
    else:
        if BaseToMutate=='A':
            NewBase = random.choice(['T','C'])
        elif BaseToMutate=='T':
            NewBase = random.choice(['A','C'])
        elif BaseToMutate=='C':
            NewBase = random.choice(['A','T'])

    NewBaseList = InputBaseList
    NewBaseList[rand_index] = NewBase
    return NewBaseList

def MutateBig(InputBaseList):
    # Input List = [W1, W2, W3, N1, N2, N3, N4, N5, N6, N7]
    rand_index = random.randint(0, 9)
    BaseToMutate = InputBaseList[rand_index]

    if rand_index < 3:
        if BaseToMutate=='A':
            NewBase = 'T'
        elif BaseToMutate=='T':
            NewBase = 'A'
    else:
        if BaseToMutate=='A':
            NewBase = random.choice(['T','C'])
        elif BaseToMutate=='T':
            NewBase = random.choice(['A','C'])
        elif BaseToMutate=='C':
            NewBase = random.choice(['A','T'])

    NewBaseList = InputBaseList
    NewBaseList[rand_index] = NewBase
    return NewBaseList

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

################################################################################
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# Sequence manipulation and mutation functions
# 
################################################################################

################################################################################
# 
# Table output functions
# \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
################################################################################

def Write_Header(filename):
    tempfile = open(filename, 'w')
    tempfile.write('Design #' + ',' +
                'W1W2' + ',' +
                'N1N2N3N4N5N6' + ',' +
                'Top spurious' + ',' +
                'Toe spurious' + ',' +
                'BM score'  + ',' +
                'Top+Toe' + ',' +
                'SPC Verboten' + ',' +
                'SPC Spurious' + ',' +
                'SPC Weighted' + ',' +
                + '\n')
    tempfile.close()

def Write_Design(Log_file_name, designcounter, Sub, TopSpuriousRel, ToeSpuriousRel, BMRel):
    tempfile = open(Log_file_name, 'a')
    tempfile.write(str(designcounter) + ',' +
                Sub['W1'] + Sub['W2'] + ',' +
                Sub['N1'] + Sub['N2'] + Sub['N3'] + Sub['N4'] + Sub['N5'] + Sub['N6'] + ',' +
                str(TopSpuriousRel) + ',' +
                str(ToeSpuriousRel) + ',' +
                str(BMRel)  + ',' +
                str(TopSpuriousRel+ToeSpuriousRel)
                + '\n')
    tempfile.close()

def HP_Write_Header(filename):
     tempfile = open(filename, 'w')
     tempfile.write('#' + ',' +
                'X1-5' + ',' +
                'Top' + ',' +
                'Max' + ',' +
                'Toe' + ',' +
                'Max' + ',' +
                'BM' + ',' +
                'Max'  + ',' +
                'SS min' + ',' +
                'SS avg' + ',' +
                '\n')
     tempfile.close()

def HP_Write_Design(Log_file_name, designcounter, Sub, Top, MaxTop, Toe, MaxToe, BM, MaxBM, SSmin, SSavg, SS_toe, SS_toeavg, ted_min, ted_avg, ted_std):
    tempfile = open(Log_file_name, 'a')
    tempfile.write(str(designcounter) + ',' + 
                str(Top) + ',' +
                str(MaxTop) + ',' +   
                str(Toe) + ',' +
                str(MaxToe) + ',' +
                str(BM)  + ',' +
                str(MaxBM)  + ',' +
                str(SSmin)  + ',' +
                str(SSavg)  + ',' +
                str(SS_toe)  + ',' +
                str(SS_toeavg)  + ',' +
                str(ted_min)  + ',' +
                str(ted_avg)  + ',' +
                str(ted_std)  + ',' +
                '\n')
    tempfile.close()

################################################################################
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# Table output functions
# 
################################################################################

################################################################################
#
# Set global peppername lists
# \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
################################################################################

def get_seq_lists(basename, gates, strands):
    # Read in sequences
    seqsfile = basename + '.seqs'
    mfefile = basename + '.mfe'
    pep_sequences, pep_strands = Read_Finished(seqsfile)
    domains_list = pep_sequences.keys()
    mfe_dict, cmplx_dict = read_design(mfefile)
    seq_dict = mfe_dict.copy()
    seq_dict.update(pep_sequences)
    seq_dict.update(pep_strands)
    
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
            if th_c not in NotToInteract:
                NotToInteract[th_c] = gate.get_noninteracting_peppernames(th)
            else:
                NotToInteract[th_c].extend(gate.get_noninteracting_peppernames(th))
            if th_c not in BaseStrandlist:
                BaseStrandlist.append(th_c)
        TopStrandlist.extend(gate.get_top_strands())
        complex_names.extend(gate.get_complexes())
        TopStranddict.update(gate.get_top_strand_dict())
    
    allpepperlist = list()
    for l in [TopStrandlist, complex_names, seq_dict.keys(), BMlist]:
        allpepperlist.extend(l)
    
    for l in NotToInteract.values():
        allpepperlist.extend(l)
    
    allpepperlist.extend([ s[:-1] for s in BaseStrandlist] )
    
    seq_dict = make_pepper_seq_dict(allpepperlist, seq_dict, update=True)
    return (TopStrandlist, complex_names, BaseStrandlist, TopStranddict, BMlist,
            NotToInteract, seq_dict, cmplx_dict, domains_list)

################################################################################
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# Set global peppername lists
# 
################################################################################

################################################################################
# 
# Score wrapper functions
# \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
################################################################################

def EvalCurrent(basename, gates, strands, compile_params=(7, 15, 2),
                header=True, testname=None):
    #import tolds_utils as tu
    # basename
    # Plan:
    #   * read MFE and SEQS file, combine?
    #   * generate pepper-name lists from species and gate objects
    #   * use these pepper-lists to make a name:sequence dictionary
    #   * pass these into necessary scoring functions...or make global?
    
    if not testname:
        testname = basename
    
    (TopStrandlist, complex_names, BaseStrandlist, TopStranddict, BMlist,
     NotToInteract, seq_dict, cmplx_dict, domains_list) = \
        get_seq_lists(testname, gates, strands)
    
    print 'Start WSI computation'
    ssm_scores = Spurious_Weighted_Score(basename, domains_list, seq_dict, 
                                         compile_params=compile_params)
    ssm_names = ['WSI-Intra', 'WSI-Inter', \
                 'WSI-Intra-1', 'WSI-Inter-1', \
                 'Verboten', 'WSI']
    print ''
    
    # Score cross-strand spurious interactions
    print 'Start Cross-Strand spurious interactions computation'
    css_scores  = NUPACK_Eval(seq_dict, TopStrandlist, BaseStrandlist, \
        NotToInteract, ComplexSize = 2, T = 25.0, material = 'dna',\
         shouldclean=1, quiet=1)
    css_names = ['TSI avg', 'TSI max', \
                 'TO avg', 'TO max']
    print ''
    
    print 'Start tube ensemble defect computation'
    ted_scores = NUPACK_Eval_tube_defect(seq_dict, cmplx_dict, complex_names,\
        prefix='tube_ensemble')
    ted_names = ['Max Bad Nucleotide %', 'Max Defect Component', 
                 'Mean Bad Nucleotide %']
    print ''
    
    # Retrieve toeholds for BM score calculation
    th_strs = [ th[:-1] for th in BaseStrandlist ]
    toeholds = [seq_dict[i] for i in th_strs]
    th_names = ['TH energy avg', 'TH energy range']
    
    # Score weighted spurious interactions
    print 'Start BM score computation'
    bm_scores = BM_Eval(seq_dict, BMlist, toeholds)
    bm_names = ['BM Score', 'Largest Match']
    print ''
    
    # Score intra-strand spurious interactions and toehold availability
    print 'Start Single-Strand spurious score computation'
    ss_scores = SS_Eval(seq_dict, TopStranddict, T = 25.0, material = 'dna')
    ss_names = ['SSU Min', 'SSU Avg', 'SSTU Min', 'SSTU Avg']
    print ''
    
    score_list = [css_scores, bm_scores, ss_scores, ted_scores, 
                  ssm_scores]
    names_list = [css_names, bm_names, ss_names, ted_names, 
                  ssm_names]
    scores = [ elem for sub in score_list for elem in sub]
    names  = [ elem for sub in names_list for elem in sub]
    
    if header:
        output = (scores, names)
    else:
        output = scores
    return output

def NUPACK_Eval_tube_defect(mfe_seqs, ideal_structs, complex_names, \
                            ComplexSize=3, T=25.0, material='dna', \
                            shouldclean=1, quiet=1, prefix='ted_calc'):
    # This function takes in the design's sequences and intended interaction 
    # structures and returns the tube ensemble defect, the concentration of incorr
    # ectly base-paired nucleotides in a tube. I think this is best prepared and t
    # hought of not as a single tube, but as if each intermediate were prepared in
    #  separate tubes and the concentration is summed over each of these tubes.
    
    # Loop through the complexes in all systems to calculate their contribution 
    # to the tube ensemble defect (ted) , the running sum of such contributions
    # Set up empty lists
    ted_vec = np.empty(len(complex_names))
    nam_list = list()
    bp_vec = np.empty(len(complex_names))
    # Set parameters
    target_conc = 1e-06
    ted_max = ComplexSize * target_conc
    max_iden = 'None'
    counter = 0
    params = [ComplexSize, T, material, quiet, prefix]
    
    prog = MyProgress(len(complex_names))
    for cmpx in complex_names:
        # Grab complex sequences and structure
        cmpx_name = cmpx
        seq_list = mfe_seqs[cmpx_name].split('+')
        seq = mfe_seqs[cmpx_name].replace('+', '')
        nseq = len(seq_list)
        struct = ideal_structs[cmpx_name]
        bp = len(re.findall('[()]', struct))
        # Retrieve the estimated complex concentration
        est_conc = NUPACK_Cmpx_Conc(seq_list, params)
        cmpx_defect = NUPACK_Cmpx_Defect(seq_list, struct, params)
        ted = cmpx_defect * min(est_conc, target_conc) + \
               len(seq) * max(target_conc - est_conc, 0)
        ted_vec[counter] = ted
        nam_list.append(cmpx_name)
        bp_vec[counter] = bp
        counter += 1
        prog.inc()
    
    # Convert TED to % bad nucleotides out of total nucleotides
    bad_nuc_vec = ted_vec / (target_conc * bp_vec + ted_vec)
    bad_nuc_max = bad_nuc_vec.max()
    bn_max_name = nam_list[np.where(bad_nuc_vec == bad_nuc_max)[0]]
    #return [Bad Nucleotide % max, max complex name, mean bad nuc]
    return [bad_nuc_max, bn_max_name, bad_nuc_vec.mean()]

def NUPACK_Eval(seq_dict, TopStrandlist, BaseStrandlist, NotToInteract,\
                ComplexSize = 2, T = 25.0, material = 'dna', \
                shouldclean=1, quiet=1):
    numstrands = len(TopStrandlist)
    
    TopSpuriousPairwise = zeros([numstrands, numstrands]);
    
    print 'Calculating Top Strand Pairwise interactions'
    countmax = (numstrands**2 + numstrands)/2
    prog = MyProgress(countmax)
    for i in range(numstrands):
        for j in range(i,numstrands):
            intij = NUPACKIntScore(TopStrandlist[i], TopStrandlist[j], 
                                   seq_dict, ComplexSize, T, material, quiet)
            TopSpuriousPairwise[i, j] = intij
            TopSpuriousPairwise[j, i] = intij
            prog.inc()
    
    numbase = len(BaseStrandlist)
    
    BaseSpurious = zeros([numbase, 1]);
    
    print 'Calculating Toehold occupation'
    countmax = sum(map(len, NotToInteract.values()))
    prog = MyProgress(countmax)
    for i in range(numbase):
        thisstrand = BaseStrandlist[i]
        notinteract = NotToInteract[thisstrand]
        nonintnum = len(notinteract)
        for j in range(nonintnum):
            intij = NUPACKIntScore(thisstrand, notinteract[j], seq_dict, 
                                   ComplexSize, T, material, quiet)
            BaseSpurious[i] = BaseSpurious[i] + intij
            prog.inc()
    
    TSI_vec = TopSpuriousPairwise.sum(0)
    return [TSI_vec.mean(), TSI_vec.max(), BaseSpurious.mean(), 
            BaseSpurious.max()]
 
    
def SS_Eval(seq_dict, TopStranddict, T = 25.0, material = 'dna'):
    TopStrandlist = TopStranddict.keys()
    numstrands = len(TopStrandlist)
    
    MinProbs = []
    avg_Unpaired = []
    
    MinProbs_toe = []
    avg_Unpaired_toe = []
    
    prog = MyProgress(numstrands)
    for i in range(numstrands):
            strand = TopStrandlist[i]
            toe_regions = TopStranddict[strand]
            [min_Unpaired, sum_Unpaired, min_Unpaired_toe, sum_Unpaired_toe,\
             NumBases] = \
                NUPACKSSScore(strand, seq_dict, T, material, toe_regions)
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
    BM_score_weights = [0,0,0,0,0,1,2,3,4,5,6,7,8,9,10,11]
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

            BM_score = BM_score + BM_score_weights[min(maxmatch, 15)]
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
                            clean=True):
    # This function calculates Niranjan's weighted spurious interaction score, for i
    # nteractions between k and km nucleotides long. It takes the following steps:
    #   * Read in sequences from mfe file
    #   * Create fixed file for all sequences except clamps
    #   * Recompile circuit
    #   * Run spuriousSSM negative design and generate scores
    #   * Grab spuriousSSM scores and calculate NSIH
      
    ###### Playin' around for spurious C scores
    import sys
    import subprocess
    import re
    import os
    from tempfile import mkstemp
    import numpy as np
    
    from sloth import call_compiler, call_design
    
    # Command parameters
    ssm_params = "bored=%s tmax=%s spurious_range=%s" % (bored, tmax, spurious_range)
    
    fid, fixed_file = mkstemp(suffix='.fixed')
    os.close(fid)
    fid, compiled_file = mkstemp(suffix='.pil')
    os.close(fid)
    fid, save_file = mkstemp(suffix='.save')
    os.close(fid)
    fid, out_file = mkstemp(suffix='.mfe')
    os.close(fid)
    fid, spurious_output = mkstemp(suffix='.txt')
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
                    fixed_file=fixed_file)
    
    # Make constraint files
    fid, design_tmp = mkstemp()
    os.close(fid)
    call_design(basename, infilename=compiled_file, outfilename=out_file, 
                just_files=True, tempname=design_tmp)
    stname = design_tmp+'.st'
    wcname = design_tmp+'.wc'
    eqname = design_tmp+'.eq'
    
    # Commented code generates MFE file
    spur_exc = 'spuriousSSM score=automatic template=%s wc=%s eq=%s %s> %s'
    command = spur_exc % (stname, wcname, eqname, ssm_params, spurious_output)
    print "Command : {}".format(command)
    os.system(command)
    #subprocess.check_call(spur_exc)
    
    w_vec = np.power(5, np.arange(spurious_range) + 1)
    
    f = open(stname)
    st = f.readline()[:-1]
    f.close()
    num_strands = len(re.findall('\s+', st))
    
    # re string for decimal/integer scraping
    num = "\d+\.?\d*"
    
    # Read in entire file, search for isntances of "spurious1"
    spc_text = open(spurious_output).read()
    vec_str = spc_text[spc_text.rfind('spurious1'):].split('\n')[2]
    vec = np.array([ np.float(x) for x in re.findall(num, vec_str)])
    score_vec = vec * w_vec
    mis_intra_score = score_vec.sum() / num_strands
    
    vec_str = spc_text[spc_text.rfind('spurious1'):].split('\n')[4]
    vec = np.array([ np.float(x) for x in re.findall(num, vec_str)])
    score_vec = vec * w_vec
    mis_inter_score = score_vec.sum() / num_strands
    
    vec_str = spc_text[spc_text.rfind('spurious('):].split('\n')[2]
    vec = np.array([ np.float(x) for x in re.findall(num, vec_str)])
    score_vec = vec * w_vec
    spc_intra_score = score_vec.sum() / num_strands
    
    vec_str = spc_text[spc_text.rfind('spurious('):].split('\n')[4]
    vec = np.array([ np.float(x) for x in re.findall(num, vec_str)])
    score_vec = vec * w_vec
    spc_inter_score = score_vec.sum() / num_strands
    
    vec_str = spc_text[spc_text.rfind('** score_verboten'):]
    verboten_score = np.float(re.findall(num, vec_str)[0]) / num_strands
    
    vec_str = spc_text[spc_text.rfind('-weighted score = '):]
    wsi_score = np.float(re.findall(num, vec_str)[-1]) / num_strands
    
    if clean:
        for f in [fixed_file, compiled_file, save_file, out_file, 
                  stname, wcname, eqname, spurious_output]:
            os.remove(f)
    
    return [spc_intra_score, spc_inter_score, \
            mis_intra_score, mis_inter_score, \
            verboten_score, wsi_score]

################################################################################
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# Score wrapper functions
# 
################################################################################

################################################################################
# 
# Nupack-calling scoring functions
# \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
################################################################################

def NUPACK_Cmpx_Conc(seqs, params=[3, 25, 'dna', 1, 'ted_calc']):
    # This function calls the NUPACK methods 'complexes' and 'concentrations' to 
    # calculate the expected concentration of the desired complex. Note that the
    # output of 'concentrations' will contain the expected concentrations of all
    # possible complexes, but this function returns only that of the desired com
    # plex, which is always 1 1 1.
    # Retrieve parameters
    from tempfile import mkstemp
    import os
    
    ComplexSize, T, material, quiet, prefix = params
    # Setup input files to nupack commands
    intsc = 0;
    fid, prefix = mkstemp()
    os.close(fid)
    ofile = prefix + '.out'
    cfile = prefix + '.con'
    ifile = prefix + '.in'
    efile = prefix + '.eq'
    file_list = [ofile, cfile, ifile, efile, prefix]
    
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
          '-material %s -ordered -pairs -mfe -dangles some -sodium 0.05 '+\
          '-magnesium 0.0125 %s > %s'
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
    
    # Count the number of commented lines
    ctr = 0;
    while (lines[ctr][0] == '%'):
        ctr = ctr + 1
    
    # Scrape the concentrations from the output file
    for count in range(ctr, len(lines)):
        lineData = whiteSpaceSearch.split(lines[count])
        if lineData[1:len(seqs)+2] == (1 + len(seqs) ) * ['1'] :
            return float(lineData[-2])
    
    for f in file_list:
        if os.path.isfile(f):
            os.remove(f)
    
    return 0


def NUPACK_Cmpx_Defect(seqs, struct, params=[3, 25, 'dna', 1, 'ted_calc']):
    '''
    This function calls the NUPACK methods 'defects' and reports the average num
    ber of basepairs that do not match their intended state.
    '''
    from tempfile import mkstemp
    import os
    # Retrieve parameters
    ComplexSize, T, material, quiet, prefix = params
    # Setup input files to nupack commands
    intsc = 0;
    fid, prefix = mkstemp()
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
    
    f.write('1 2 3\n')
    f.write(struct)
    f.close()
    
    # Run NUPACK's defect function
    cmd = nupackpath + 'defect -T %.1f '+\
          '-material %s -dangles some -multi %s > %s'
    cmd = cmd % (T,material,prefix,ofile) 
    os.system(cmd)
    
    # Read the equilibrium concentrations 
    f = open(ofile)
    lines = f.readlines()
    f.close()
    
    # Clean
    for f in file_list:
        if os.path.isfile(f):
            os.remove(f)
    
    # Count the number of commented lines
    ctr = 0
    while (lines[ctr][0] == '%'):
        ctr = ctr + 1
    return float(lines[ctr])


def NUPACKIntScore(str1, str2, seq_dict, ComplexSize, T, material, quiet):
    # Formerly IntScore
    from tempfile import mkstemp
    import os
    
    intsc = 0;
    fid, fname = mkstemp()
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
          '-material %s -ordered -pairs -mfe -dangles some -sodium 0.05 '+\
          '-magnesium 0.0125 %s > %s' 
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
    
    # Clean
    for f in file_list:
        if os.path.isfile(f):
            os.remove(f)

    for count in range(ctr, len(lines)):
        lineData = whiteSpaceSearch.split(lines[count])
        rawsc = float(lineData[-2])
        intsc = intsc + 100*rawsc/(2*1e-6)
    
    return intsc

def NUPACKSSScore(str1, seq_dict, T, material, toe_region=[None]):
    # This function calls NUPACK's complexes and concentrations functions, which
    # allows one to calculate the unpaired probabilities of each nucleotide across
    # the most common secondary structures.
    from tempfile import mkstemp
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
    fid, fname = mkstemp()
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
          '-material %s -ordered -pairs -mfe -dangles some -sodium 0.05 '+\
          '-magnesium 0.0125 %s > %s' 
    cmd = cmd % (T,material,fname,ofile)
    os.system(cmd)
    cmd = nupackpath + 'concentrations  '+\
          '-ordered -pairs -cutoff 0.000001 -sort 0 -quiet ' + fname
    os.system(cmd)
    
    # Read in nucleotide pairing information
    f = open(pfile)
    lines = f.readlines()
    f.close()
    
    # Clean
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
        


################################################################################
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# Nupack-calling scoring functions
# 
################################################################################
