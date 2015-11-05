from __future__ import division
import sys
import numpy as np

# from DSDClasses import *
# from LeaklessClasses import *

def call_compiler(basename, 
                    args = (7, 15, 2), 
                    outputname=None, 
                    savename=None, 
                    fixed_file=None, 
                    synth=True, 
                    includes=
                      ['/home/jmp/Dropbox/CRNtoDNASoftware/crn_compiler']):
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

def write_toehold_file(toehold_file, strands, toeholds):
    """ Writes the fixed file for the given strands and toeholds
    
    Args:
        toehold_file: File where toehold constraints are written
        strands: List of SignalStrand objects
        toeholds: List of toehold two-big tuples
    
    Returns:
        Nothing
    """
    line = 'sequence {} = {} # species {}\n'
    f = open(toehold_file, 'w')
    for strand, ths in zip(strands, toeholds):
        for i in range(len(ths)):
            constraint = line.format(strand.th(i), ths[i].upper(), 
                                     strand.get_names()[-1])
            f.write(constraint)
    f.close()

def set_signal_instances_constraints(basename, strands, mod):
    """ Wrapper for writing instance equivalence 'reactions' to the sys file
    
    Args:
        basename: The default name of filetypes to be produced and accessed.
        strands: List of SignalStrand objects
        mod: module containing scheme variables and classes
    
    Returns:
        i : reaction index
    """
    i = 0
    sysfile = basename + '.sys'
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
    basename = sysfile[:-4]
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

def write_compiler_files(basename, gates, strands, toeholds, mod):
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
    sysfile = basename + '.sys'
    fixedfile = basename + '.fixed'
    write_sys_header(sysfile, mod)
    
    # Write sys file
    f = open(sysfile, 'a')
    for gate in gates:
        f.write(gate.get_reaction_line())
    f.close()
    
    # append signal strand constraints to the sys file and write toehold fixed
    # file
    set_signal_instances_constraints(basename, strands, mod)
    write_toehold_file(fixedfile, strands, toeholds)

def generate_scheme(basename, th_params=dict(), design_params=(7, 15, 2), \
                    fixedfile=None, mod=None):
    """ Produce PIL file for a CRN system
    
    A scheme consists toeholds and a PIL file. In the standard workflow, a
    single scheme will produce a set of sequences conforming to the parameters
    and constraints of the scheme. Generating a scheme requires producing the
    constraints and system file describing a DSD implementation of a CRN. This
    function, making use of other methods, reads in a text file describing an 
    abstract CRN, determines the strands and gates necessary to implement the
    CRN, assigns toeholds to signal strands, and produces input files for
    downstream sequence design processes.
        
    Args:
        basename: Default name for files accessed and written
        th_params: Dictionary of toehold design parameters ( dict() )
            thold_l: Nt in a toehold (7)
            thold_e: Target deltaG in kCal/Mol (7.7)
            e_dev: Allowable standard deviation in kCal/mole (0.5)
            m_spurious: Maximum spurious dG as fraction of thold_e (0.4)
            e_module: Thermodynamics used by stickydesign (energyfuncs_james)
        design_params: A tuple of parameters to the system file ( (7, 15, 2) )
        fixedfile: Filename for the sequence constraints (None)
        mod: module containing scheme variables and classes (DSDClasses)
    Returns:
        gates: A list of gate objects
        strands: A list of strand objects
        toeholds: A list of two-tuples containing paired toeholds
        th_scores: Average and range of toehold binding energies
    """
    if mod is None:
        import DSDClasses as mod
    
    crn_file = basename + ".crn"
    systemfile = basename + ".sys"
    pilfile = basename + ".pil"
    
    if fixedfile==None:
        fixedfile = basename + ".fixed"
    else:
        # Expect user to have produced an appropriate fixed file
        pass
    
    reactions, species = read_crn(crn_file)
    
    output = mod.process_rxns(reactions, species, design_params)
    (gates, strands) = output
    
    n_species = len(strands)
    signal_names = [ s.get_names()[-1] for s in strands ]
    labs = [' first', ' second']
    labels = [ spec + lab for spec in signal_names for lab in labs ]
    toeholds, th_scores = toehold_wrapper(n_species*mod.n_th, th_params, labels)
    toehold_sets = [ toeholds[i:i+mod.n_th] for i in range(0,mod.n_th*n_species,mod.n_th)]
    write_compiler_files(basename, gates, strands, toehold_sets, mod) 
    try:
        call_compiler(basename, args=design_params, fixed_file=fixedfile)
    except Exception, e:
        print e
    
    return (gates, strands, toeholds, th_scores)

def generate_seqs(basename, 
                  design_params=(7, 15, 2), 
                  outname=None,
                  extra_pars=""):
    """ Produce sequences file for scheme
    
    This function finds a file -basename.pil- and produces sequences using
    spuriousSSM negative design.
    
    Args:
        basename: Default name for files accessed and written
        design_params: A tuple of parameters to the system file ( (7, 15, 2) )
        outname: Optional name for files produced as a result of this call
        extra_pars: Options sent to spurious designer. ('')
    Returns:
        Nothing
    """
    pilfile = basename + ".pil"
    mfefile = basename + ".mfe"
    seqfile = basename + ".seqs"
    strandsfile = basename + "_strands.txt"
    
    if outname:
        mfefile = outname + '.mfe'
        seqfile = outname + ".seqs"
        strandsfile = outname + "_strands.txt"
    
    # Now do the sequence makin' 
    call_design(basename, pilfile, mfefile, verbose=True, 
                extra_pars=extra_pars)
    
    # "Finish" the sequence generation
    call_finish(basename, designname=mfefile, seqsname=seqfile, run_kin=False)

def run_designer(basename='small', 
                 reps=1, 
                 th_params={"thold_l":7, "thold_e":7.7, "e_dev":1, \
                            "m_spurious":0.5, "e_module":'energyfuncs_james'},
                 design_params=(7, 15, 2),
                 mod_str=None,
                 extra_pars=""):
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
    Returns:
        Nothing, but writes many basename + extension files, such as:
            system file (.sys)
            sequences (.seqs)
            scores (_scores.csv)
    """
    if mod_str is None:
        mod_str = 'DSDClasses'
    
    exec('import {} as mod'.format(mod_str))
    import tdm
    import numpy as np
    fixedfile = basename + ".fixed"
    systemfile = basename + ".sys"
    pilfile = basename + ".pil"
    
    (gates, strands, toeholds, th_scores) = \
        generate_scheme(basename, th_params, design_params, mod=mod)
    
    if reps => 1:
        trialnames = [ basename + str(i) for i in range(reps) ]
        for testname in trialnames:
            try:
                generate_seqs(basename, design_params, testname, extra_pars=extra_pars)
            except Exception as e:
                print 'Error!'
                print e
                return (gates, strands)
        
        scoreslist = []
        for testname in trialnames:
            try:
                scores, score_names = tdm.EvalCurrent(basename, gates, strands, 
                                                  testname=testname, 
                                                  compile_params=design_params)
            except Exception as e:
                print 'Error!'
                print e
                return (gates, strands)
            scoreslist.append(scores)
        
        score_names = ['Set Index', 'Toehold Avg dG', 'Range of toehold dG\'s'] + score_names
        for i in range(len(scoreslist)):
            l = scoreslist[i]
            scoreslist[i] = [i, th_scores[0], th_scores[1]] + l
        
        f = open(basename+'_scores.csv', 'w')
        f.write(','.join(score_names))
        f.write('\n')
        f.writelines( [ ','.join(map(str, l)) + '\n' for l in scoreslist ])
        f.close()
        outstring = [ score_names[i+1] + ':{}  '.format(scoreslist[i]) for i in range(len(scoreslist))]
        print outstring
    return (gates, strands)

def run_secretary_problem(basename='small', 
                         reps=2, 
                         th_params={"thold_l":7, "thold_e":7.7, "e_dev":1, \
                                    "m_spurious":0.5, "e_module":'energyfuncs_james'},
                         design_params=(7, 15, 2),
                         mod_str=None,
                         extra_pars="bored=10"):
    # Setup
    if mod_str is None:
        mod_str = 'DSDClasses'
    
    exec('import {} as mod'.format(mod_str))
    import tdm
    import numpy as np
    fixedfile = basename + ".fixed"
    systemfile = basename + ".sys"
    pilfile = basename + ".pil"
    
    (gates, strands, toeholds, th_scores) = \
        generate_scheme(basename, th_params, design_params, mod=mod)
    
    # Begin
    if reps > 1:
        # Generate burn-through candidates
        trialnames = [ basename + str(i) for i in range(reps) ]
        for testname in trialnames:
            try:
                generate_seqs(basename, design_params, testname, extra_pars=extra_pars)
            except Exception as e:
                print 'Error!'
                print e
                return (gates, strands)
        
        # Score burn-through candidates
        scoreslist = []
        score_comp = np.array([1,1,1,1,1,1,-1,-1,-1,-1,1,1,1,1,1,1,1,1])
        for testname in trialnames:
            try:
                scores, score_names = tdm.EvalCurrent(basename, gates, strands, 
                                                  testname=testname, 
                                                  compile_params=design_params)
                score_vec = np.array(scores)
                score_vec = score_vec[np.concatenate((np.arange(11), np.arange(12, len(score_vec))))]
                scores = scores * score_comp
            except Exception as e:
                print 'Error!'
                print e
                return (gates, strands)
            scoreslist.append(scores)
        
        # Determine benchmark
        score_arr = np.array(scoreslist)
        score_arr.sort(0)
        score_best = np.mean(score_arr, 0)
        print "best score"
        print score_best
        
        # Generate candidates until benchmark is met
        while np.any(score - score_best < 0):
            print "testing new talent"
            testname = "selection"
            generate_seqs(basename, design_params, testname, extra_pars=extra_pars)
            score, score_names = tdm.EvalCurrent(basename, gates, strands, 
                                              testname=testname, 
                                              compile_params=design_params)
            score = np.array(score)
            score = score[np.concatenate((np.arange(11), np.arange(12, len(score))))]
            scores = score * score_comp
            
        print "winning score"
        print score
        scoreslist.append(score)
        score_names = ['Set Index', 'Toehold Avg dG', 'Range of toehold dG\'s'] + score_names
        for i in range(len(scoreslist)):
            l = scoreslist[i]
            scoreslist[i] = [i, th_scores[0], th_scores[1]] + l
        
        f = open(basename+'_scores.csv', 'w')
        f.write(','.join(score_names))
        f.write('\n')
        f.writelines( [ ','.join(map(str, l)) + '\n' for l in scoreslist ])
        f.close()
        outstring = [ score_names[i+1] + ':{}  '.format(scoreslist[i]) for i in range(len(scoreslist))]
        print outstring
    return 

def score_fixed(fixedfile, 
                 basename='small', 
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
    crn_file = basename + ".crn"
    sysfile = basename + ".sys"
    pilfile = basename + ".pil"
    extra_pars = ""
    
    reactions, species = read_crn(crn_file)
    output = mod.process_rxns(reactions, species, design_params)
    (gates, strands) = output
    
    # Write header immediately
    write_sys_header(sysfile, mod)
    
    # Write sys file
    f = open(sysfile, 'a')
    for gate in gates:
        f.write(gate.get_reaction_line())
    f.close()
    
    # append signal strand constraints to the sys file and write toehold fixed
    # file
    set_signal_instances_constraints(basename, strands, mod)
    try:
        call_compiler(basename, args=design_params, fixed_file=fixedfile)
    except Exception, e:
        print e
    
    generate_seqs(basename, design_params, extra_pars=extra_pars)
    scores, score_names = tdm.EvalCurrent(basename, gates, strands, 
                                      compile_params=design_params)
    return (scores, score_names)

if __name__ == "__main__":
    mod_str = 'DSDClasses'
    #mod_str = 'LeaklessClasses'
    basename = 'oscillator'
    reps = 1
    th_params = {"thold_l":7, "thold_e":7.7, "e_dev":1, \
                 "m_spurious":0.5, "e_module":'energyfuncs_james'}
    t = 7
    bm = 15
    c = 2
    design_params = (t,bm,c)
    #design_params = (7,)
    gates, strands = run_designer(basename, reps, th_params, design_params, mod_str, 
                                    extra_pars="bored=10")
