from __future__ import division
import sys
import numpy as np
import os.path

def call_compiler(basename, 
                    args = (7, 15, 2), 
                    outputname=None, 
                    savename=None, 
                    fixed_file=None, 
                    synth=True, 
                    includes=
                      [os.path.dirname(__file__)]):
    """ Generates a PIL file from a .sys and fixed files.
    
    Args:
        basename: The default name of filetypes to be produced and accessed.
        args: A tuple of system arguments. 
        outputname: the PIL file produced by the compiler. def: basename.pil
        savename: the save file produced by the compiler. def: basename.save
        fixed_file: filename specifying sequence constraints. def: none
        synth: Boolean, whether or not to produce an output. Deprecated.
        includes: path to folder holding component files referenced by .sys
    Returns:
        Nothing
    """
    from PepperCompiler.compiler import compiler
    if outputname is None:
        outputname = '{}.pil'.format(basename)
    if savename is None:
        savename = '{}.save'.format(basename)
    if os.path.dirname(basename) is not '':
        includes = includes + [os.path.dirname(basename)]
        basename = os.path.basename(basename)
    compiler(basename, args, outputname, savename, fixed_file, synth, includes)

def call_design(basename, 
                infilename=None, 
                outfilename=None, 
                cleanup=True, 
                verbose=False, 
                reuse=False, 
                just_files=False, 
                struct_orient=False, 
                old_output=False, 
                tempname=None, 
                extra_pars="", 
                findmfe=True, 
                spuriousbinary="spuriousSSM"):
    """ Generates an MFE file from a .pil file.
    
    Args:
        basename: The default name of filetypes to be produced and accessed.
        infilename: PIL file read by the designer. (basename.pil)
        outfilename: MFE file read by the designer. (basename.mfe)
        cleanup: Boolean, Delete st wc eq files
        verbose: Boolean, verbose designer output
        reuse: Boolean, if available and appropriate, use wc st eq files
        just_files: Boolean, only produce wc st eq files
        struct_orient: Boolean, list structures in MFE file
        old_output: Deprecated
        tempname: Optional temporary name for wc st eq files
        extra_pars: Options sent to spurious designer. ('')
        findmfe: Use DNAfold to do something mysterious. (True)
        spuriousbinary: Compiled C++ for negative design. (spuriousSSM)
    Returns:
        Nothing
    """
    from PepperCompiler.design import spurious_design as spd
    if not infilename:
        infilename = '{}.pil'.format(basename)
    if not outfilename:
        outfilename = '{}.mfe'.format(basename)
    spd.design(basename, infilename, outfilename, cleanup, verbose, reuse, 
              just_files, struct_orient, old_output, tempname, extra_pars, 
              findmfe, spuriousbinary)

def call_finish(basename,
                savename=None,
                designname=None,
                seqsname=None,
                strandsname=None,
                run_kin=False,
                cleanup=True,
                trials=24,
                time=1e6,
                temp=25.0,
                conc=1,
                spurious=False,
                spurious_time=10.0):
    """ Generates a .seqs file from an .mfe file.
    
    Args:
        basename: The default name of filetypes to be produced and accessed.
        savename: File storing process states. (basename.save)
        designname: MFE file read for sequences (basename.mfe)
        seqsname: Output file containing all sequences (basename.seqs)
        strandsname: Output file containing all strand sequences
        run_kin: Run spurious kinetic tests on sequences (False)
        cleanup: Delete temporary files (True)
        trials: Number of kinetics trials to run (24)
        time: Simulation seconds (1000000)
        temp: Degrees celsius for simulations (25)
        conc: Concentration of strands in simulations in uM (1)
        spurious: Run pairwise kinetics tests (false)
        spurious_time: Simulation time (s) for spurious kinetics (10.0)
    Returns:
        Nothing
    """
    from PepperCompiler import finish
    if not savename:
        savename = '{}.save'.format(basename)
    if not designname:
        designname = '{}.mfe'.format(basename)
    if not seqsname:
        seqsname = '{}.seqs'.format(basename)
    finish.finish(savename, designname, seqsname, strandsname, run_kin, 
                  cleanup, trials, time, temp, conc, spurious, spurious_time)

def read_crn(in_file):
    """ Interprets a CRN text file.
    
    Maintains a list of reactions and species. For each reaction, the rate 
    constant, reactants, products, and stoichiometric coefficients are read.
    The interpreter relies on a lot of regex to extract tokens from the lines
    of texts. 
    
    Tokens are as follows:
        * Alphanumeric are species
        * Numeric before alphanumeric are stoichiometric identifiers
        * Arrow (->) separates reactants from products
        * Parentheses enclose the rate constant, which python should be able 
          to evaluate
    Args:
        in_file : String that this the file name to be read
    
    Returns:
        crn_info : A tuple containing lists 'reactions' and 'species'. 
                   'reactions' contains a list of dicionaries keyed by
                   'reactants', 'products', 'rate', 'stoich_r', and 
                   'stoich_p'
    """
    import re 
    fid = open(in_file, 'r')
    lines = list()
    for line in fid:
        lines.append(line[:-1])
    
    fid.close()
    
    # Read through the file and extract reaction specifications
    rxn_tup = list()
    spe_ind_dic = dict()
    species_list = list()
    update_ind = 0
    num_pattern = re.compile(r"^[0-9./]+|^[0-9.]+e-?[0-9.]+")
    spe_pattern = re.compile(r"\w+")
    for line in lines:
        # eat empty lines
        if line == '':
            continue
        
        full_line = line[:]
        # Remove whitespace
        line = re.sub(r'\s','',line)
        
        # Check for reaction rate
        rate_match = re.search(r'\(.*\)', line)
        if rate_match:
            rate_str = rate_match.group(0)[1:-1]
            rate = float(eval(rate_str))
            line = re.sub(r'\(' + rate_str + '\)', '', line)
        else:
            rate = 1
        
        # Split into reactants and products
        rxn_eq = re.split('->',line)
        try:
            assert len(rxn_eq)==2
        except AssertionError, e:
            print 'Check reaction arrow use.'
            raise
            
        stoich = 1
        reactants = list()
        products = list()
        stoich_r = list()
        stoich_p = list()
        for subline, lhs_flag in zip(rxn_eq, (True, False)):
            tokens = re.split(r'[ +]+', subline)
            for token in tokens:
                # Attempt to scrape for coefficients
                coeff_string = num_pattern.search(token)
                if coeff_string:
                    coeff_string = num_pattern.search(token).group(0)
                    try:
                        stoich = eval(coeff_string)
                    except SyntaxError, e:
                        print 'improper coefficient syntax!'
                        raise
                    token = token.replace(coeff_string, '')
                else:
                    stoich = 1
                
                # Remaining token string should be species identifier
                if spe_pattern.match(token):
                    if token not in spe_ind_dic:
                        spe_ind_dic[token] = update_ind
                        update_ind += 1
                        species_list.append(token)
                    
                    spe = spe_ind_dic[token]
                    
                    if lhs_flag:
                        reactants.append(token)
                        stoich_r.append(stoich)
                    else:
                        products.append(token)
                        stoich_p.append(stoich)
                    
                    stoich = 1
                    
        rxn_tup.append({"reactants":reactants,
                        "products":products, 
                        "stoich_r":stoich_r, 
                        "stoich_p":stoich_p, 
                        "rate":rate})
    
    return (rxn_tup, species_list)

def list_to_str(in_specs):
    # Write a single string composed of the strings in the input list separated
    # by plus signs and spaces
    # early return for empty case
    if len(in_specs) == 0:
        return ''
    # Make a deep copy to allow the pops
    specs = in_specs[:]
    spec_str = specs.pop()
    if specs:
        for spec in specs:
            spec_str = spec_str + " + " + spec
    return spec_str

def write_signal_equals(sysfile, species, i, mod):
    """ Writes a line to sysfile setting two signal identities equal
    
    The set_signal_equals component takes two species (writen as A -> B) and
    sets the toes and identity branch migration domains euqal. This is used to
    set the different produce instances of signals equal to one another yet
    allow them different history domains
    
    Args:
        sysfile: The system file to write the reaction line to
        species: Species instances to be set equal
        i: reaction index
        mod: module containing scheme variables and classes
    
    Returns:
        i: increased reaction index
    """
    comp = mod.set_equals
    f = open(sysfile, 'a')
    lines = []
    base_line = '{}eq{} = {}{}: {} -> {}\n'
    if len(species) > 2:
        spec_a = species.pop()
        for spec_b in species:
            line = base_line.format(mod.rxn_strings[0], str(i), comp, 
                                    mod.rxn_strings[1], spec_a, spec_b)
            lines.append(line)
            i = i + 1
    elif len(species) == 2:
        spec_a = species[0]
        spec_b = species[1]
        line = base_line.format(mod.rxn_strings[0], str(i), comp, 
                                mod.rxn_strings[1], spec_a, spec_b)
        lines.append(line)
        i = i + 1
    else:
        f.close()
        return i
    f.writelines(lines)
    f.close()
    return i

def write_toehold_file(toehold_file, strands, toeholds, n_th):
    """ Writes the fixed file for the given strands and toeholds
    
    Args:
        toehold_file: File where toehold constraints are written
        strands: List of SignalStrand objects
        toeholds: List of toehold two-big tuples
    
    Returns:
        Nothing
    """
    from numpy import mod as modulo
    from numpy import floor 
    line = 'sequence {} = {} # species {}\n'
    f = open(toehold_file, 'w')
    if type(toeholds[0]) is str:
        for i, th in zip(range(len(toeholds)), toeholds):
            strand = strands[ int(floor((i - modulo(i, n_th)) / n_th)) ]
            constraint = line.format(strand.th(modulo(i, n_th)), toeholds[i].upper(), 
                                     strand.get_names()[-1])
            f.write(constraint)
    else:
        for strand, ths in zip(strands, toeholds):
            for i in range(len(ths)):
                constraint = line.format(strand.th(i), ths[i].upper(), 
                                         strand.get_names()[-1])
                f.write(constraint)
    f.close()

def set_signal_instances_constraints(sysfile, strands, mod):
    """ Wrapper for writing instance equivalence 'reactions' to the sys file
    
    Args:
        basename: The default name of filetypes to be produced and accessed.
        strands: List of SignalStrand objects
        mod: module containing scheme variables and classes
    
    Returns:
        i : reaction index
    """
    i = 0
    species = []
    for strand in strands:
        instances = strand.get_names()
        if len(instances) > 2: # species is a product more than once
            i = write_signal_equals(sysfile, instances, i, mod)
    return i

def write_sys_header(sysfile, mod):
    """ Writes header for the CRN system file
    
    Writes the system declaration line and import statements. No need to be
    intelligent about the import statements used, so this is a pretty dumb
    function
    
    Args:
        sysfile: System file passed to PepperCompiler
        mod: module containing scheme variables and classes
    
    Returns:
        Nothing
    """
    basename = os.path.basename(sysfile[:-4])
    f = open(sysfile, 'w')
    f.write("declare system " + basename + mod.param_string + " -> \n")
    f.write("\n")
    # Comps is defined in Classes file
    for comp in mod.comps:
        f.write("import {0}\n".format(comp))
    f.write("\n")
    f.close()

def toehold_wrapper(n_ths, paramdict=dict(), labels=None):
    """ Wrapper generating toeholds, calls gen_th
    
    Args:
        n_spec: Number of species requiring two toeholds each
        paramdict: Dictionary of parameters used in producing toeholds
            thold_l: Nt in a toehold (7)
            thold_e: Target deltaG in kCal/Mol (7.7)
            e_dev: Allowable standard deviation in kCal/mole (0.5)
            m_spurious: Maximum spurious dG as fraction of thold_e (0.4)
            e_module: Thermodynamics used by stickydesign (energyfuncs_james)
    Returns:
        ths: Toeholds, listed as tupled-pairs
        th_score: average toehold dG and the range of dG's
    """
    from gen_th import get_toeholds
    # Grab parameters from the dictionary or set defaults
    # Toehold length (basepairs)
    if 'thold_l' in paramdict:
        thold_l = paramdict['thold_l']
    else:
        thold_l = int(7)
    # Target toehold energy
    if 'thold_e' in paramdict:
        thold_e = paramdict['thold_e']
    else:
        thold_e = 7.7
    # Allowable energy deviation in kcal/mol
    if 'e_dev' in paramdict:
        e_dev = paramdict['e_dev']
    else:
        e_dev = 0.5
    # max spurious interaction as a fraction of target energy
    if 'm_spurious' in paramdict:
        m_spurious = paramdict['m_spurious']
    else:
        m_spurious = 0.4
    # Energetics class to use in thermodynamics calculations
    if 'e_module' in paramdict:
        e_module = paramdict['e_module']
    else:
        e_module = 'energyfuncs_james'
    thold_l = int(thold_l)
    ths, th_score =  \
        get_toeholds(n_ths, thold_l, thold_e, e_dev, m_spurious, e_module)
    return (ths, th_score)

def write_compiler_files(basename, 
                         gates, 
                         strands,
                         mod,
                         sysfile=None,
                         inst_constraints=True):
    """ Write input files to the Pepper compiler

    This function iterates through each reaction recording the gates and
    strands required to model it using DNA. These strand and gate objects
    hold references to specific strands, complexes, domains, and subsequences
    of the DSD system that help construct the inputs to the Pepper compiler 
    and spuriousSSM. In addition, they allow the scoring functions to access 
    the specific nucleotide sequences they require.
    
    Args:
        basename: Default name for files accessed and written
        gates: A list of gate objects
        strands: A list of strand objects
        toeholds: A list of two-big tuples holding toehold strings
        mod: module containing scheme variables and classes
    Returns:
        Nothing
    """
    
    # Write header immediately
    if sysfile is None:
        sysfile = basename + '.sys'
    write_sys_header(sysfile, mod)
    
    # Write sys file
    f = open(sysfile, 'a')
    for gate in gates:
        f.write(gate.get_reaction_line())
    f.close()
    
    # append signal strand constraints to the sys file and write toehold fixed
    # file
    if inst_constraints:
        set_signal_instances_constraints(sysfile, strands, mod)

def generate_scheme(basename,  
                    design_params=(7, 15, 2), 
                    mod=None,
                    crn_file=None,
                    systemfile=None,
                    inst_constraints=True):
    """ Produce SYS file describing a CRN
    
    A scheme consists a .sys file and lists of gate and strand objects. Gate 
    and strand objects tell the scoring modules, write_compiler_files, and 
    write_toehold_files which names in the .PIL file refer to the DNA sequence 
    domains these functions need to access.  This function, making use of 
    other methods, reads in a text file describing an abstract CRN, determines 
    the strands and gates necessary to implement the CRN, and writes a .sys 
    for a DNA approximation of that reaction network.
        
    Args:
        basename: Default name for files accessed and written
        design_params: A tuple of parameters to the system file ( (7, 15, 2) )
        mod: module containing scheme variables and classes (DSDClasses)
    Returns:
        gates: A list of gate objects
        strands: A list of strand objects
    """
    if mod is None:
        import DSDClasses as mod
    
    if crn_file is None:
        crn_file = basename + ".crn"
    if systemfile is None:
        systemfile = basename + ".sys"
    
    reactions, species = read_crn(crn_file)
    
    output = mod.process_rxns(reactions, species, design_params)
    (gates, strands) = output
    
    write_compiler_files(basename, gates, strands, mod, systemfile, inst_constraints) 
    return (gates, strands)

def generate_seqs(basename, 
                  gates, 
                  strands, 
                  design_params=(7, 15, 2), 
                  n_th=2, 
                  th_params={"thold_l":7, "thold_e":7.7, "e_dev":1, \
                             "m_spurious":0.5, "e_module":'energyfuncs_james'},
                  outname=None,
                  extra_pars="",
                  systemfile=None,
                  pilfile=None,
                  mfefile=None,
                  seqfile=None,
                  fixedfile=None,
                  savefile=None,
                  strandsfile=None):
    """ Produce sequences for a scheme
    
    This function accepts a base file name, a list of gate objects, a list of 
    strand objects, and toehold parameters and calls StickyDesign, the 
    PepperCompiler, and SpuriousSSM to generate a DNA sequence.     
    
    Args:
        basename: Default name for files accessed and written
        gates: A list of Gate objects
        strands: A list of SignalStrand objects
        design_params: A tuple of parameters to the system file ( (7, 15, 2) )
        n_th: How many toeholds to generate per signal strand
        th_params: Dictionary of toehold design parameters 
            thold_l: Nt in a toehold (7)
            thold_e: Target deltaG in kCal/Mol (7.7)
            e_dev: Allowable standard deviation in kCal/mole (0.5)
            m_spurious: Maximum spurious dG as fraction of thold_e (0.4)
            e_module: Thermodynamics used by stickydesign (energyfuncs_james)
        outname: Optional name for files produced as a result of this call
        extra_pars: Options sent to spurious designer. ('')
    Returns:
        Nothing
    """
    
    # Prepare filenames
    if systemfile is None:
        systemfile = basename + ".sys"
    if pilfile is None:
        pilfile = basename + ".pil"
    if savefile is None:
        savefile = basename + ".save"
    if mfefile is None:
        if outname:
            mfefile = outname + '.mfe'
        else:
            mfefile = basename + ".mfe"
    if seqfile is None:
        if outname:
            seqfile = outname + ".seqs"
        else:
            seqfile = basename + ".seqs"
    if strandsfile is None:
        if outname:
            strandsfile = outname + "_strands.txt"
        else:
            strandsfile = basename + "_strands.txt"
    if fixedfile is None:
        fixedfile = basename + ".fixed"
    
    # Make toeholds
    n_species = len(strands)
    signal_names = [ s.get_names()[-1] for s in strands ]
    labs = [' first', ' second']
    labels = [ spec + lab for spec in signal_names for lab in labs ]
    toeholds, th_scores = toehold_wrapper(n_species*n_th, th_params, labels)
    toehold_sets = [ toeholds[i:(i+n_th)] for i in range(0,n_th*n_species,n_th)]
    # Write the fixed file for the toehold sequences and compile the sys file to PIL
    write_toehold_file(fixedfile, strands, toeholds, n_th)
    try:
        call_compiler(basename, args=design_params, fixed_file=fixedfile, 
                      outputname=pilfile, savename=savefile)
    except Exception, e:
        print e
    
    # Now do the sequence makin' 
    call_design(basename, pilfile, mfefile, verbose=True, 
                extra_pars=extra_pars, cleanup=False)
    
    # "Finish" the sequence generation
    call_finish(basename, savename=savefile, designname=mfefile, \
                seqsname=seqfile, strandsname=strandsfile, run_kin=False)
    return (toeholds, th_scores)

def selection(scores):
    columns = zip(*scores)
    ranks = []
    fractions = []
    percents = []
    for col in columns:
        # print col
        if 'Index' in col[0] or 'Defect' in col[0] or 'Toehold Avg' in col[0] or 'Range of toehold' in col[0]:
            continue
        if 'SSU' in col[0] or 'SSTU' in col[0]:    # for these scores, higher is better
            col = [-float(x) for x in col[1:]]
        else:
            col = [float(x) for x in col[1:]]
        array = np.array(col)
        temp = array.argsort()
        colranks = np.empty(len(array), int)
        colranks[temp] = np.arange(len(array))
        ranks.append(colranks)                     # low rank is better
        fractions.append((array - array.min())/abs(array.min() + (array.min()==0) ))
        percents.append((array - array.min())/(array.max() - array.min()))
    temp=ranks
    ranks=zip(*temp)
    temp=fractions
    fractions=zip(*temp)
    temp=percents
    percents=zip(*temp)
    
    print "\nRank array:"
    print "                         ",
    for title in scores[0]:
        if 'Index' in title or 'Defect' in title or 'Toehold Avg' in title or 'Range of toehold' in title:
            continue
        print "{:>6s}".format(title[0:6]),
    print
    i=0
    for r in ranks:
        print "design {:2d}: {:6d} = sum [".format(i,sum(r)),
        for v in r:
            print "{:6d}".format(v),
        print "]"
        i=i+1
    
    print "\nFractional excess array:"
    print "                         ",
    for title in scores[0]:
        if 'Index' in title or 'Defect' in title or 'Toehold Avg' in title or 'Range of toehold' in title:
            continue
        print "{:>6s}".format(title[0:6]),
    print
    i=0
    for f in fractions:
        print "design {:2d}: {:6.2f} = sum [".format(i,sum(f)),
        for v in f:
            print "{:6.2f}".format(v),
        print "]"
        i=i+1
    
    print "\nPercent badness (best to worst) array:"
    print "                         ",
    for title in scores[0]:
        if 'Index' in title or 'Defect' in title or 'Toehold Avg' in title or 'Range of toehold' in title:
            continue
        print "{:>6s}".format(title[0:6]),
    print
    i=0
    for p in percents:
        print "design {:2d}: {:6.2f} = sum [".format(i,100*sum(p)),
        for v in p:
            print "{:6.2f}".format(100*v),
        print "]"
        i=i+1
    
    print " "
    worst_rank = 0
    while 1:
        ok_seqs = [i for i in range(len(ranks)) if max(ranks[i])<=worst_rank]
        if len(ok_seqs)==0:
            worst_rank=worst_rank+1
            continue
        else:
            break
    
    # scores used:
    # TSI avg, TSI max, TO avg, TO max, BM, Largest Match, SSU Min, SSU Avg, SSTU Min, SSTU Avg, Max Bad Nt %,  Mean Bad Nt %, WSI-Intra, WSI-Inter, WSI-Intra-1, WSI-Inter-1, Verboten, WSI
    weights = [5,   20,     10,     30,  2,             3,      30,      10,       50,       20,           10,              5,         6,         4,           5,           3,        2,  8] 
        
    print "Indices of sequences with best worst rank of " + str(worst_rank) + ": " + str(ok_seqs)
    print "  Sum of all ranks, for these sequences:      " + str([sum(ranks[i]) for i in ok_seqs])
    print "  Sum of weighted ranks, for these sequences: " + str([sum(np.array(ranks[i])*weights/100.0) for i in ok_seqs])
    print "  Sum of fractional excess over best score:   " + str([sum(fractions[i]) for i in ok_seqs])
    print "  Sum of weighted fractional excess:          " + str([sum(np.array(fractions[i])*weights/100.0) for i in ok_seqs])
    print "  Sum of percent badness scores:              " + str([100*sum(percents[i]) for i in ok_seqs])
    print "  Sum of weighted percent badness scores:     " + str([sum(np.array(percents[i])*weights) for i in ok_seqs])
    temp = [sum(r) for r in ranks]
    print "Best sum-of-ranks:                   {:6.2f} by [{:d}]      and the worst: {:6.2f} by [{:d}]".format( min(temp), np.argmin(temp), max(temp), np.argmax(temp) )
    winner = np.argmin(temp)
    temp = [sum(np.array(r)*weights/100.0) for r in ranks]
    print "Best sum-of-weighted-ranks:          {:6.2f} by [{:d}]      and the worst: {:6.2f} by [{:d}]".format( min(temp), np.argmin(temp), max(temp), np.argmax(temp) )
    temp = [sum(f) for f in fractions]
    print "Best fractional excess sum:          {:6.2f} by [{:d}]      and the worst: {:6.2f} by [{:d}]".format( min(temp), np.argmin(temp), max(temp), np.argmax(temp) )
    temp = [sum(np.array(f)*weights/100.0) for f in fractions]
    print "Best weighted fractional excess sum: {:6.2f} by [{:d}]      and the worst: {:6.2f} by [{:d}]".format( min(temp), np.argmin(temp), max(temp), np.argmax(temp) )
    temp = [100*sum(p) for p in percents]
    print "Best percent badness sum:            {:6.2f} by [{:d}]      and the worst: {:6.2f} by [{:d}]".format( min(temp), np.argmin(temp), max(temp), np.argmax(temp) )
    temp = [sum(np.array(p)*weights) for p in percents]
    print "Best weighted percent badness sum:   {:6.2f} by [{:d}]      and the worst: {:6.2f} by [{:d}]".format( min(temp), np.argmin(temp), max(temp), np.argmax(temp) )
    print
    return winner

def selection_wrapper(scores, reportfile = 'score_report.txt'):
    import sys
    stdout = sys.stdout
    try:
        sys.stdout = open(reportfile, 'w')
        winner = selection(scores)
    except Exception, e:
        sys.stdout.close()
        sys.stdout = stdout
        print e
        print sys.exc_info()[0]
        raise
        return 'bad'
    else:
        sys.stdout.close()
        sys.stdout = stdout
    return winner

def run_designer(basename='small', 
                 reps=1, 
                 th_params={"thold_l":7, "thold_e":7.7, "e_dev":1, \
                            "m_spurious":0.5, "e_module":'energyfuncs_james'},
                 design_params=(7, 15, 2),
                 mod=None,
                 extra_pars="",
                 quick=False):
    """ Generate and score sequences
    
    This function links together all component of the compiler pipeline. It 
    generates a system and PIL file, then runs negative sequence design 
    software multiple times to generate multiple sequences. The scoring wrappper
    accepts the gates and strands lists of objects that allow the scoring
    functions to access and compute on the domains and strands from each
    sequence set.
    
    Args:
        basename: Default name for files accessed and written
        reps: Number of sequence sets to be generated and scored (2)
        th_params: Dictionary of toehold design parameters 
            thold_l: Nt in a toehold (7)
            thold_e: Target deltaG in kCal/Mol (7.7)
            e_dev: Allowable standard deviation in kCal/mole (0.5)
            m_spurious: Maximum spurious dG as fraction of thold_e (0.4)
            e_module: Thermodynamics used by stickydesign (energyfuncs_james)
        design_params: A tuple of parameters to the system file ( (7, 15, 2) )
        mod_str: A string of the name of the .py file containing the scheme 
                 definitions
        extra_pars: Options sent to spurious designer. ('')
        quick: Make random scores instead of computing heursitics. Skips time
               consuming computations for debugging purposes. (False)
    Returns:
        Nothing, but writes many basename + extension files, such as:
            system file (.sys)
            sequences (.seqs)
            scores (_scores.csv)
    """
    if mod is None:
        import DSDClasses as mod
    
    import tdm
    fixedfile = basename + ".fixed"
    systemfile = basename + ".sys"
    pilfile = basename + ".pil"
    
    (gates, strands) = \
        generate_scheme(basename, design_params, mod)
    
    if reps >= 1:
        scoreslist = []
        for i in range(reps):
            testname = basename + str(i) + '.txt'
            try:
                toeholds, th_scores = generate_seqs(basename, gates, strands, design_params, \
                              mod.n_th, th_params, strandsfile=testname, extra_pars=extra_pars)
                scores, score_names = tdm.EvalCurrent(basename, gates, strands, 
                                                  testname=testname, 
                                                  compile_params=design_params,
                                                  quick=quick)
                scores = [i] + [s for s in th_scores] + [s for s in scores]
                scoreslist.append(scores)
            except Exception as e:
                print 'Error!'
                print e
                return (gates, strands)
        score_names = ['Set Index', 'Toehold Avg dG', 'Range of toehold dG\'s'] + score_names
        scores = [score_names] + scoreslist
        winner = selection_wrapper(scores)
        
        f = open(basename+'_scores.csv', 'w')
        f.write(','.join(score_names))
        f.write('\n')
        f.writelines( [ ','.join(map(str, l)) + '\n' for l in scoreslist ])
        f.close()
    return (gates, strands, winner)

def score_fixed(fixed_file, 
                 basename='small', 
                 crn_file=None, 
                 sys_file=None, 
                 pil_file=None, 
                 save_file=None, 
                 mfe_file=None, 
                 seq_file=None,
                 design_params=(7, 15, 2),
                 mod_str=None):
    """ Score a sequence set
    
    This function takes in a fixed file, crn file, and reaction scheme specification
    and outputs the heuristic scores of the set. 
    
    Args:
        fixedfile: Filename pointing to the seqeunce set
        basename: Default name for files accessed and written
        design_params: A tuple of parameters to the system file ( (7, 15, 2) )
        mod_str: A string of the name of the .py file containing the scheme 
                 definitions
    Returns:
        scores: A list containing the scores generated by EvalCurrent
        score_names: A list of strings describing the scores
    """
    if mod_str is None:
        mod_str = 'DSDClasses'
    
    exec('import {} as mod'.format(mod_str))
    import tdm
    if crn_file is None:
        crn_file = basename + '.crn'
    if sys_file is None:
        sys_file  = basename + '.sys'
    if pil_file is None:
        pil_file = basename + '.pil'
    if save_file is None:
        save_file = basename + '.save'
    if mfe_file is None:
        mfe_file = basename + '.mfe'
    if seq_file is None:
        seq_file = basename + '.seq'
    extra_pars = ""
    
    gates, strands = generate_scheme(basename,  
                                   design_params, 
                                   crn_file=crn_file, 
                                   systemfile=sys_file)
    try:
        call_compiler(basename, args=design_params, fixed_file=fixed_file, 
                      outputname=pil_file, savename=save_file)
    except Exception, e:
        print e
    
    # Now do the sequence makin' 
    call_design(basename, pil_file, mfe_file, verbose=True, 
                extra_pars=extra_pars, cleanup=False)
    # "Finish" the sequence generation
    call_finish(basename, savename=save_file, designname=mfe_file, \
                seqsname=seq_file, run_kin=False)
    scores, score_names = tdm.EvalCurrent(basename, gates, strands, 
                                      compile_params=design_params, 
                                      seqs_file=seq_file, mfe_file=mfe_file)
    return (scores, score_names)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-b", "--basename", help='Basename.[small]', type=str)
    parser.add_argument("-l", "--length", help='Toehold length.[7]', type=int)
    parser.add_argument("-e", "--energy", help='Target toehold binding energy'+
                        'in kcal/mol.[7.7]', type=float)
    parser.add_argument("-d", "--deviation", help='Toehold binding energy st'+\
                        'andard deviation limit in kcal/mol.[0.5]', type=float)
    parser.add_argument("-M", "--maxspurious", help='Maximum spurious interac'+
                        'tion energy as a multiple of target binding energy.['+
                        '0.4]', type=float)
    parser.add_argument("-g", '--energetics', help='The name of the energetics'+
                        ' module to be used for toehold generation.[energyfunc'+
                        's_james]', type=str)
    parser.add_argument("-p", '--systemparams', help='A string of integers that'+
                        ' are parameters to the sys file compilation[Finds in module]', type=str)
    parser.add_argument("-n", '--candidates', help='Number of candidate sequences'+
                        ' to generate[1]', type=int)
    parser.add_argument("-m", '--module', help='Module describing the strand'+
                        ' displacement architecture[DSDClasses]', type=str)
    parser.add_argument("-x", '--extrapars', help='Parameters sent to SpuriousSSM[]', type=str)
    parser.add_argument("-q", '--quick', action='store_true', 
                        help='Make random numbers instead of computing heuristics to save time[False]')
    args = parser.parse_args()
    ############## Interpret arguments
    if args.basename:
        basename = args.basename
    else:
        basename = 'small'
    # Set toehold length in base pairs
    if args.length:
        thold_l = int(args.length)
    else:
        thold_l = int(7)
    
    # Set target toehold binding energy, kcal/mol
    if args.energy:
        thold_e = args.energy
    else:
        thold_e = 7.7
    
    # Set allowable toehold binding energy deviation, kcal/mol
    if args.deviation:
        e_dev = args.deviation
    else:
        e_dev = 0.5
    
    # Set allowable spurious interactions, fraction of target energy.
    # The spurious interactions are calculated with standard SantaLucia parameters
    # and sticky-end contexts, not toehold contexts
    if args.maxspurious:
        m_spurious = args.maxspurious
    else:
        m_spurious = 0.4
    
    # Import desired energetics module as ef
    if args.energetics:
        energetics = args.energetics
    else:
        energetics = 'energyfuncs_james'
    
    # Set user-specified spurious interaction heatmap image file name
    if args.candidates:
        reps = args.candidates
    else:
        reps = 1
    
    if args.module:
        exec('import {} as mod'.format(args.module))
    else:
        import DSDClasses as mod
    
    # Set user-specified spurious interaction heatmap image file name
    if args.systemparams:
        design_params = [ int(i) for i in args.systemparams.split(' ')]
    else:
        try:
            design_params = mod.default_params
        except Exception, e:
            print e
            print 'Cannot guess default .sys parameters without a module'
            raise
    
    if args.extrapars:
        extra_pars = args.extrapars
    else:
        extra_pars = None
    
    th_params={"thold_l":thold_l, "thold_e":thold_e, "e_dev":e_dev, \
               "m_spurious":m_spurious, "e_module":energetics},
    
    gates, strands, winner = run_designer(basename, reps, th_params, design_params, mod, 
                                    extra_pars=extra_pars, quick=args.quick)
    print 'Winning sequence set is index {}'.format(winner)
