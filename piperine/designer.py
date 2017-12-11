#!/usr/bin/env python
from __future__ import division, print_function

import sys
import os
import os.path
import importlib
import pkg_resources

import numpy as np

from . import Srinivas2017 as default_translation_scheme

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

default_energyfuncs = default_translation_scheme.energetics.energyfuncs()
default_translation = default_translation_scheme.translation
default_design_params = default_translation_scheme.translation.default_params
small_crn = pkg_resources.resource_filename('piperine', "data/small.crn")
data_dir = os.path.dirname(small_crn)

def call_compiler(basename,
                  args=default_design_params,
                  outputname=None,
                  savename=None,
                  fixed_file=None,
                  includes=None):
    """ Generates a PIL file from a .sys. This is a wrapper for a peppercompiler function.

    Args:
        basename: The default name of filetypes to be produced and accessed.
        args: A tuple of system arguments.
        outputname: The PIL file produced by the compiler. (Default: <basename>.pil)
        savename: The save file produced by the compiler. def: (Default: <basename>.save)
        fixed_file: Filename specifying sequence constraints. (Default: No fixed file)
        includes: Path to folder holding component files referenced by .sys.
                  (Default: piperine.data)
    Returns:
        Nothing
    """
    from peppercompiler.compiler import compiler
    if outputname is None:
        outputname = '{}.pil'.format(basename)
    if savename is None:
        savename = '{}.save'.format(basename)
    if includes is None:
        includes = []
    includes.append(data_dir)
    compiler(basename, args, outputname, savename, fixed_file, True, includes)

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
                findmfe=False,
                spuriousbinary="spuriousSSM"):
    """ Generates an MFE file from a .pil file. This is a wrapper for a peppercompiler function.

    Args:
        basename: The default name of filetypes to be produced and accessed.
        infilename: PIL file read by the designer. (basename.pil)
        outfilename: MFE file read by the designer. (basename.mfe)
        cleanup: Boolean, Delete st wc eq files (True)
        verbose: Boolean, verbose designer output (True)
        reuse: Boolean, if available and appropriate, use wc st eq files (False)
        just_files: Boolean, only produce wc st eq files (False)
        struct_orient: Boolean, list structures in MFE file (False)
        old_output: Deprecated (False)
        tempname: Optional temporary name for wc st eq files (None)
        extra_pars: Options sent to spurious designer. ('')
        findmfe: Use DNAfold to do something mysterious. (True)
        spuriousbinary: Compiled C++ for negative design. (spuriousSSM)
    Returns:
        Nothing
    """
    from peppercompiler.design.spurious_design import design
    if not infilename:
        infilename = '{}.pil'.format(basename)
    if not outfilename:
        outfilename = '{}.mfe'.format(basename)
    design(basename, infilename, outfilename, cleanup, verbose, reuse,
           just_files, struct_orient, old_output, tempname, extra_pars,
           findmfe, spuriousbinary)
    if not os.path.isfile(outfilename):
        raise RuntimeError('Expected MFE not created, expect SSM failure')

def call_finish(basename,
                savename=None,
                designname=None,
                seqname=None,
                strandsname=None,
                run_kin=False,
                cleanup=True,
                trials=24,
                time=1e6,
                temp=25.0,
                conc=1,
                spurious=False,
                spurious_time=10.0):
    """ Generates a .seqs file from an .mfe file. This is a wrapper for a peppercompiler function.

    Args:
        basename: The default name of all files produced and accessed.
        savename: File storing process states. (Default: basename.save)
        designname: MFE file, read for sequences (Default: basename.mfe)
        seqname: Output file containing all sequences (Default: basename.seqs)
        strandsname: Output file containing all strand sequences (Default: None)
        run_kin: Run spurious kinetic tests on sequences (Default: False)
        cleanup: Delete temporary files (Default: True)
        trials: Number of kinetics trials to run (Default: 24)
        time: Simulation seconds (Default: 1000000)
        temp: Degrees celsius for simulations (Default: 25)
        conc: Concentration of strands in simulations in uM (Default: 1)
        spurious: Run pairwise kinetics tests (Default: False)
        spurious_time: Simulation time (s) for spurious kinetics (Default: 10.0)
    Returns:
        Nothing
    """
    from peppercompiler.finish import finish
    if savename is None:
        savename = '{}.save'.format(basename)
    if designname is None:
        designname = '{}.mfe'.format(basename)
    if seqname is None:
        seqname = '{}.seqs'.format(basename)

    finish(savename, designname, seqname, strandsname, run_kin,
           cleanup, trials, time, temp, conc, spurious, spurious_time)

def parse_parameter_line(line, translation=default_translation):
    """ This interprets compilation or DNA domain parameters found in the CRN file.

    This function is called whenever the get_parameters_from_crn_file function
    encounters an equals sign "=" and, if the parameter term is recognized, the
    associated value is recorded and assigned internally to the matched parameter.

    Parameter terms are either those related to toehold generation, defining
    translation scheme package, sequence designer The following are parameters that can
    be assigned through this method. Terms in parentheticals are the terms this
    function detects.

        * Any parameters required by components (.comp) files. The translation scheme
          package provides the terms to look for.
        * Toehold parameters
            * (toehold_energy) Target binding energy
            * (toehold_deviation) Maximum standard deviation of toehold binding energies
            * (toehold_spurious) Maximum spurious binding energy
        * (translation scheme) Python package defining translation and energyfuncs
        * (spurious_design_parameters) A string that is appended to the
          spurious design call. Look at the spurious design documentation for info.
        * (n) Number of sets to generate

    Other terms will be accepted if they appear in the list of parameter terms
    defined in the translation scheme. Terms not found in these places will
    be ignored. Terms are detected in line order, so higher line number is
    higher priority.

    Args:
        line : String from the CRN file defining a compilation parameter
        translation : Translation module (Default: Srinivas2017.translation)

    Returns:
        param_dict : A dictionary listing the parameter definitions found in the CRN file
    """
    import re
    # Accepted terms
    terms = ["toehold_energy",
             "toehold_deviation",
             "toehold_spurious",
             "toehold_length",
             "spurious_design_parameters",
             "translation_scheme",
             "n"]

    converters = [lambda x: np.float(x),
                  lambda x: np.float(x),
                  lambda x: np.float(x),
                  lambda x: np.int(x),
                  lambda x: x,
                  lambda x: x,
                  lambda x: np.int(x)]

    # Design parameters are always integers
    if  translation is not None:
        for param_term in default_translation.param_terms:
            terms.append(param_term)
            converters.append(lambda x: np.int(x))

    # Remove whitespace
    line = re.sub(r'\s', '', line)
    rhs, lhs = line.split('=')

    param_dict = {}
    for i, term in enumerate(terms):
        if rhs == term:
            param_dict.update({term:converters[i](lhs)})
    return param_dict

def get_parameters_form_crn_file(crn_file, translation=default_translation):
    '''
    Wraapper function for parse_parameter_line.

    Args:
        line : String from the CRN file defining a compilation parameter
        translation : Translation module (Default: Srinivas2017.translation)

    Returns:
        param_dict : A dictionary listing the parameter definitions found in the CRN file
    '''
    fid = open(crn_file, 'r')
    lines = list()
    for line in fid:
        lines.append(line[:-1])

    fid.close()

    # Read through the file and extract reaction specifications
    parameters = {}
    for line in lines:
        if '=' in line:
            parameters.update(parse_parameter_line(line, translation))

    return parameters

def read_crn(in_file):
    """ Interprets a CRN from a text file.

    Maintains a list of reactions and species. For each reaction, the rate
    constant, reactants, products, and stoichiometric coefficients are read.
    The interpreter relies on a lot of regex to extract tokens from the lines
    of texts.

    Tokens are as follows:
        * Alphanumeric are species
        * Numeric before alphanumeric are stoichiometric identifiers
        * Arrow (->) separates reactants from products
        * Lines with an equals signs (=) are skipped
        * Parentheses enclose the rate constant, which python should be able
          to evaluate
    Args:
        in_file : String of the file name holding the CRN specification

    Returns:
        crn_info : A tuple containing lists called 'reactions' and 'species'.
                   'reactions' contains a list of dicionaries keyed by
                   'reactants', 'products', 'rate', 'stoich_r', and
                   'stoich_p'. 'species' list holds strings representing
                   signal species names.
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
    parameters = {}
    for line in lines:
        # Skip empty lines
        if line == '':
            continue

        # Skip lines that define design parameters
        if '=' in line:
            continue

        full_line = line[:]
        # Remove whitespace
        line = re.sub(r'\s', '', line)

        # Check for reaction rate
        rate_match = re.search(r'\(.*\)', line)
        if rate_match:
            rate_str = rate_match.group(0)[1:-1]
            rate = float(eval(rate_str))
            line = re.sub(r'\(' + rate_str + r'\)', '', line)
        else:
            rate = 1

        # Split into reactants and products
        rxn_eq = re.split('->', line)
        try:
            assert len(rxn_eq) == 2
        except AssertionError as e:
            print('Check reaction arrow use.')
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
                    except SyntaxError as e:
                        print('improper coefficient syntax!')
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
    # Build list of PIL names for toehold domains
    th_names = [th for strand in strands for th in strand.get_ths()]
    th_names = sorted(list(set(th_names)))
    th_data = [(th, ', '.join([strand.name for strand in strands if th in strand.get_ths()]))
                    for th in th_names]
    # Write file
    with open(toehold_file, 'w') as f:
        for data, seq in zip(th_data, toeholds):
            constraint = line.format(data[0], seq.upper(), data[1])
            f.write(constraint)
        f.close()

def write_sys_file(basename,
                   gates=None,
                   sys_file=None,
                   translation=default_translation):
    """ Write system file from gates list

    This function takes in filenames and a CRN specification and writes a system
    file, which outlines the overall DNA implementation.

    Args:
        basename: Default name for files accessed and written
        gates: List of gate objects
        sys_file: System file filename (basename + .sys)
        translation: module containing scheme variables and classes (Srinivas2017.translation)
    Returns:
        Nothing
    """

    # Write header immediately
    if sys_file is None:
        sys_file = basename + '.sys'

    # Clean basename if it was provided with directory prefix
    if os.path.sep in basename:
        basename = os.path.basename(basename)

    with open(sys_file, 'w') as f:
        f.write("declare system " + basename + translation.param_string + " -> \n")
        f.write("\n")
        # Comps is defined in Classes file
        for comp in translation.comps:
            f.write("import {0}\n".format(comp))
        f.write("\n")
        for rxn in gates:
            f.write(rxn.get_reaction_line())

def process_crn(basename=None,
                design_params=None,
                translation=default_translation,
                crn_file=None):
    """ Generate objects describing DNA implementation

    Gate and strand objects tell the scoring modules, write_sys_file, and
    write_toehold_files which names in the .PIL file refer to the DNA sequence
    domains these functions need to access.  This function, making use of
    other methods, reads in a text file describing an abstract CRN and determines
    the strands and gates necessary to implement the CRN.
    Args:
        basename: Default name for files accessed and written
        design_params: A tuple of parameters to the system file ( (7, 15, 2) )
        translation: module containing scheme variables and classes (Srinivas2017.translation)
        crn_file: name of the text file specifying the CRN (basename + .crn)
    Returns:
        gates: A list of gate objects
        strands: A list of strand objects
    """

    if crn_file is None:
        crn_file = basename + ".crn"

    if design_params is None:
        design_params = translation.default_params

    reactions, species = read_crn(crn_file)

    output = translation.process_rxns(reactions, species, design_params)
    (gates, strands) = output

    return (gates, strands)

def generate_scheme(basename,
                    design_params=default_design_params,
                    translation=default_translation,
                    crn_file=None,
                    system_file=None):
    """ Produce SYS file describing a CRN

    A scheme consists a .sys file and lists of gate and strand objects. Gate
    and strand objects tell the scoring modules, write_sys_file, and
    write_toehold_files which names in the .PIL file refer to the DNA sequence
    domains these functions need to access.  This function, making use of
    other methods, reads in a text file describing an abstract CRN, determines
    the strands and gates necessary to implement the CRN, and writes a .sys
    for a DNA approximation of that reaction network.

    Args:
        basename: Default name for files accessed and written
        design_params: A tuple of parameters to the system file ( (7, 15, 2) )
        translation: module containing scheme variables and classes (Srinivas2017.translation)
        crn_file: name of the text file specifying the CRN (basename + .crn)
        system_file: name of the system file (basename + .sys)
    Returns:
        gates: A list of gate objects
        strands: A list of strand objects
    """

    if system_file is None:
        system_file = basename + ".sys"

    (gates, strands) = process_crn(basename, design_params, translation, crn_file)

    write_sys_file(basename, gates, system_file, translation)
    return (gates, strands)

def generate_seqs(basename,
                  gates,
                  strands,
                  design_params=default_design_params,
                  energyfuncs=default_energyfuncs,
                  outname=None,
                  extra_pars="",
                  system_file=None,
                  pil_file=None,
                  mfe_file=None,
                  seq_file=None,
                  fixed_file=None,
                  save_file=None,
                  strands_file=None):
    """ Produce sequences for a scheme

    This function accepts a base file name, a list of gate objects, a list of
    strand objects, and toehold parameters and calls StickyDesign, the
    peppercompiler, and SpuriousSSM to generate a DNA sequence.

    Args:
        basename: Default name for files accessed and written
        gates: A list of Gate objects
        strands: A list of SignalStrand objects
        design_params: A tuple of parameters to the system file ( (7, 15, 2) )
        energyfuncs: Thermodynamics used by stickydesign
        outname: Optional name for files produced as a result of this call (basename + .mfe)
        extra_pars: Options sent to spurious designer. ("")
        system_file: Filename of the peppercompiler system file (basename + .sys)
        pil_file: Filename of the peppercompiler (basename + .pil)
        mfe_file: Filename of the peppercompiler MFE file (basename + .mfe)
        seq_file: Filename of the peppercompiler sequence file (basename + .seqs)
        fixed_file: Filename of the peppercompiler fixed file (Default: None)
        save_file: Filename of the peppercompiler save file (basename + .save)
        strands_file: Filename of the peppercompiler strands file (basename + _strands.txt)
    Returns:
        toeholds:
    """

    # Prepare filenames
    if system_file is None:
        system_file = basename + ".sys"
    if pil_file is None:
        pil_file = basename + ".pil"
    if save_file is None:
        save_file = basename + ".save"
    if mfe_file is None:
        if outname is not None:
            mfe_file = outname + '.mfe'
        else:
            mfe_file = basename + ".mfe"
    if seq_file is None:
        if outname:
            seq_file = outname + ".seqs"
        else:
            seq_file = basename + ".seqs"
    if strands_file is None:
        if outname:
            strands_file = outname + "_strands.txt"
        else:
            strands_file = basename + "_strands.txt"
    if fixed_file is None:
        fixed_file = basename + ".fixed"

    # Make toeholds
    n_species = len(strands)
    signal_names = [s.name for s in strands]
    tdomains = []
    for strand in strands:
        tdomains += strand.get_ths()
    tdomains = list(set(tdomains))
    n_toeholds = len(tdomains)

    toeholds = energyfuncs.get_toeholds(n_toeholds)

    # Write the fixed file for the toehold sequences and compile the sys file to PIL
    write_toehold_file(fixed_file, strands, toeholds)
    try:
        with Capturing() as cptr:
            call_compiler(basename, args=design_params, fixed_file=fixed_file,
                          outputname=pil_file, savename=save_file)
    except KeyError as e:
        raise(e)

    # Generate sequences
    with Capturing() as cptr:
        call_design(basename, pil_file, mfe_file, verbose=False,
                    extra_pars=extra_pars, cleanup=False)

    # "Finish" the sequence generation
    with Capturing() as cptr:
        call_finish(basename, savename=save_file, designname=mfe_file, \
                    seqname=seq_file, strandsname=strands_file, run_kin=False)
    print('Sequences written to: {}'.format(mfe_file))
    print('Strand identities written to: {}'.format(strands_file))
    return toeholds

def selection_wrapper(scores, reportfile = 'score_report.txt', optimizer='metarank'):
    '''
    Calculate relative optimality based on various rank-based statistics.

    Args:
        scores: List of list of scores from each candidate sets.
        reportfile: Filename for the output of the score report table ('score_report.txt')
        optimizer: Nonparametric measure used to select the optimum candidate ('metarank')
    Returns:
        Winning sequence set index
    '''
    from . import selectseq
    stdout = sys.stdout
    if optimizer == 'ranksum':
        selection = lambda x : selectseq.ranksum(x)
    elif optimizer == 'metarank':
        selection = lambda x : selectseq.metarank(x)
    else:
        pass
        "TD: Put error message here"
    try:
        sys.stdout = open(reportfile, 'w')
        winner = selection(scores)
    except Exception as e:
        sys.stdout.close()
        sys.stdout = stdout
        print(e)
        print(sys.exc_info()[0])
        raise
    else:
        sys.stdout.close()
        sys.stdout = stdout
    return winner

def run_designer(basename,
                 reps=4,
                 design_params=default_design_params,
                 translation=default_translation,
                 energyfuncs=default_energyfuncs,
                 optimizer='metarank',
                 extra_pars="",
                 quick=False,
                 includes=None
                ):
    """ Generate and score sequences

    This function links together all component of the compiler pipeline. It
    generates a system and PIL file, then runs negative sequence design
    software multiple times to generate multiple sequences. The scoring wrappper
    accepts the gates and strands lists of objects that allow the scoring
    functions to access and compute on the domains and strands from each
    sequence set.

    Writes many basename + extension files, such as:
        system file (.sys)
        sequences (.seqs)
        scores (_scores.csv)

    Args:
        basename: Default name for files accessed and written (small)
        reps: Number of sequence sets to be generated and scored (1)
        design_params: A tuple of parameters to the system file ( (7, 15, 2) )
        translation: Module containing scheme variables and classes (Srinivas2017.translation)
        energyfuncs: Module defining all functions related to toehold generation (Srinivas2017.energetics())
        extra_pars: Options sent to spurious designer. ('')
        quick: Make random scores instead of computing heursitics. Skips time
               consuming computations for debugging purposes. (False)
        includes: List of directories that peppercompiler looks in to find .comp files (None)
    Returns:
        gates: List of gate objects created by the translation module
        strands: List of strand objects created by the translation module
        winner: Integer defining the index of the winning candidate
        scoreslist: List of list containing scores for each candidate
    """
    # Provide extra parameters to spuriousSSM such that no minimization is performed
    if quick:
        extra_pars = "imax=-1 quiet=TRUE"

    from . import tdm
    system_file = basename + ".sys"
    pil_file = basename + ".pil"

    (gates, strands) = \
        generate_scheme(basename, design_params, translation)

    if reps >= 1:
        scoreslist = []
        for i in range(reps):
            testname = basename + str(i) + '.txt'
            seqsname = basename + str(i) + '.seqs'
            try:
                toeholds = generate_seqs(basename,
                                         gates,
                                         strands,
                                         design_params,
                                         energyfuncs=energyfuncs,
                                         strands_file=testname,
                                         seq_file=seqsname,
                                         extra_pars=extra_pars)

                scores, score_names = tdm.EvalCurrent(basename,
                                                      gates,
                                                      strands,
                                                      testname=testname,
                                                      seq_file=seqsname,
                                                      compile_params=design_params,
                                                      quick=quick,
                                                      includes=includes,
                                                      energyfuncs=energyfuncs)
                scores = [i] + scores
                scoreslist.append(scores)
            except KeyError as e:
                print('Error!')
                print(e)
                return (gates, strands, e)
        score_names = ['Set Index'] + score_names
        scores = [score_names] + scoreslist
        if reps > 2:
            winner = selection_wrapper(scores,
                       reportfile=basename+'_score_report.txt', optimizer=optimizer)
        else:
            winner = None
        with open(basename+'_scores.csv', 'w') as f:
            f.write(','.join(score_names))
            f.write('\n')
            f.writelines( [','.join(map(str, l)) + '\n' for l in scoreslist])

    return (gates, strands, winner, scoreslist)

def score_fixed(fixed_file,
                 basename=os.path.dirname(__file__)+'/small',
                 crn_file=None,
                 sys_file=None,
                 pil_file=None,
                 save_file=None,
                 mfe_file=None,
                 seq_file=None,
                 score_file=None,
                 design_params=default_design_params,
                 translation=default_translation,
                 energyfuncs=default_energyfuncs,
                 includes=None,
                 quick=False):
    """ Score a sequence set

    This function takes in a fixed file, crn file, and reaction scheme specification
    and outputs the heuristic scores of the set.

    Args:
        fixed_file: Filename pointing to the seqeunce set
        basename: Default name for files accessed and written
        crn_file: Filename of text file specifying the CRN (basename + .crn)
        sys_file: Filename of the peppercompiler system file (basename + .sys)
        pil_file: Filename of the peppercompiler PIL file (basename + .pil)
        save_file: Filename of the peppercompiler save file (basename + .save)
        mfe_file: Filename of the peppercompiler MFE file (basename + .mfe)
        seq_file: Filename of the peppercompiler output seq file (basename + .seqs)
        design_params: A tuple of parameters to the system file ( (7, 15, 2) )
        translation: Module containing scheme variables and classes (Srinivas2017.translation)
        optimizer: Method used by the selection process to choose the optimal candidate. ('metarank')
        energyfuncs: Module defining all functions related to toehold generation (Srinivas2017.energetics())
        includes: List of directories that peppercompiler looks in to find .comp files (None)
        quick: Skip time-consuming steps of minimizing sequence symetry and scoring (False)
    Returns:
        scores: A list containing the scores generated by EvalCurrent
        score_names: A list of strings describing the scores
    """
    from . import tdm
    if crn_file is None:
        crn_file = basename + '.crn'
    else:
        basename_temp = crn_file[:-4]
    if sys_file is None:
        sys_file  = basename + '.sys'
    else:
        basename_temp = sys_file[:-4]
    if pil_file is None:
        pil_file = basename + '.pil'
    else:
        basename_temp = pil_file[:-4]
    if save_file is None:
        save_file = basename + '.save'
    else:
        basename_temp = save_file[:-5]
    if mfe_file is None:
        mfe_file = basename + '.mfe'
    else:
        basename_temp = mfe_file[:-4]
    if seq_file is None:
        seq_file = basename + '.seqs'
    else:
        basename_temp = seq_file[:-4]
    if score_file is None:
        score_file = fixed_file[:-6] + '.score'
    if basename is None:
        basename = bn
    extra_pars = ""

    # Generate .sys file
    gates, strands = generate_scheme(basename,
                                   design_params,
                                   crn_file=crn_file,
                                   system_file=sys_file,
                                   translation=translation)
    # Generate PIL
    with Capturing() as output:
        call_compiler(basename, args=design_params, fixed_file=fixed_file,
                      outputname=pil_file, savename=save_file, includes=includes)

    # Generate MFE
    with Capturing() as output:
        call_design(basename, pil_file, mfe_file, verbose=False,
                    extra_pars=extra_pars, cleanup=False)

    # Generate .seqs file
    with Capturing() as output:
        call_finish(basename, savename=save_file, designname=mfe_file, \
                    seqname=seq_file, run_kin=False, cleanup=False)

    scores, score_names = tdm.EvalCurrent(basename,
                                          gates,
                                          strands,
                                          compile_params=design_params,
                                          quick=quick,
                                          includes=includes,
                                          seq_file=seq_file,
                                          mfe_file=mfe_file,
                                          energyfuncs=energyfuncs)
    with open(score_file, 'w') as f:
        f.write(','.join(score_names))
        f.write('\n')
        f.writelines([','.join(map(str, l)) + '\n' for l in [scores]])
    return (scores, score_names)

def get_design_parser():
    import argparse
    descr = '\
piperine-design is the command-line utility for compiling CRNs into DNA sequences using Piperine. From this utility, \
users can specify the domain length parameters, toehold energetics, and translation scheme. The outputs of this utility \
are the sequences of candidate DNA implementations, a table of heuristic measures of sequence quality, and a comparison \
of the candidates based on their heuristic scores. The utility will also announce the candidate that Piperine considers \
to be the best.'

    usage = "\n" + descr + "\n\n\n\
The following is the call template with short option flags. Options are shown in brackets. Capitalized terms stand in for \
arguments. Ellipsis indicate multiple arguments accepted.\n\n\
piperine-design CRNFILE [-h] [-l LENGTH] [-e ENERGY] [-d DEVIATION] [-m MAXSPURIOUS] \
[-p DESIGNPARAMS ...] [-n CANDIDATES] \
[-t TRANSLATION] [-x EXTRAPARS] [-q] \n\n\
Default execution parameters are stated in the option descriptions. The following are example executions with CRN \
file my_very_own.crn and option arguments to override the default settings. Files generated will have the same file name \
as the .crn file, but have different extensions (e.g. my_very_own.pil, my_very_own0.seqs, my_very_own_score_report.txt). \
Sequence candidates will be generated in \n\
    \n\
    Design one set of sequences for the my_very_own.crn according to the default translation scheme, Srinivas2017:\n \n\
    piperine-design my_very_own.crn\n \
\n \
\n\
    Design 10 sequence sets for my_very_own.crn:\n\n\
    piperine-design my_very_own.crn -n 10\n\
\n \
\n\
    Override default domain lengths for the default translation scheme:\n\n\
    piperine-design my_very_own.crn --designparams 5 20 2 \n\
\n \
\n\
    Design one set of sequences for my_very_own.crn according to the Chen2013 translation scheme:\n\n\
    piperine-design my_very_own.crn -t Chen2013\n\
    or\n\
    piperine-design my_very_own.crn --translation_scheme Chen2013\n\
    \n \n\
    Override default domain lengths for the Chen2013 translation scheme:\n\n\
    piperine-design my_very_own.crn --designparams 5 20 --translation_scheme Chen2013\n\n\
Parameters may also be set in the CRN file using equality statements and the following terms \
(e.g. toehold_energy=7.5):\n\
    - toehold_energy\n\
    - toehold_deviation\n\
    - toehold_spurious\n\
    - toehold_length\n\
    - spurious_design_parameters\n\
    - translation_scheme\n\
    - n\n"
    parser = argparse.ArgumentParser(usage=usage)
    parser.add_argument("crnfile",
                        help='Text file describing CRN. May also define parameters. Files read or generated by ' +
                              'this call will have the form: basename.extension ; '+
                          'E.g. basename_score_report.txt, basename.pil, when the crn file is basename.crn.',
                        type=str)

    parser.add_argument("-l", "--length",
                        help='Integer specifying toehold length in nucleotides. [default: Taken from translation scheme package]',
                        type=int)

    parser.add_argument("-e", "--energy",
                        help='Target toehold binding energy used by StickyDesign, '+
                        'in kcal/mol. [default: 7.7]',
                        type=float)

    parser.add_argument("-d", "--deviation",
                        help='Maximum standard deviation for toehold binding energies allowed by StickyDesign'+
                            ', in kcal/mol. [default: 0.5]',
                        type=float)

    parser.add_argument("-m", "--maxspurious",
                        help='This argument is passed to StickyDesign and sets the Maximum spurious interac'+
                        'tion energy as a multiple of target binding energy. [default: '+
                        '0.4]',
                        type=float)

    parser.add_argument("-p", '--designparams',
                        help='A string of integers that are parameters to the sys file compilation. [default: Finds in '+
                             'translation scheme module]',
                        type=int,
                        nargs="*")

    parser.add_argument("-n", '--candidates',
                        help='Number of candidate sequences to generate. [default: 4]',
                        type=int)

    parser.add_argument("-O", '--optimizer',
                        help='Method used by the selection process to choose the optimal candidate. Options are "ranksum" or \
                              "metarank". [default: metarank]',
                        type=str)

    parser.add_argument("-t", '--translation_scheme',
                        help='Provide a string, the name of the Python package describing the translation scheme used to convert'+
                            'the CRN to DNA strands and complexes. See piperine.Srinivas2017 for an example of such a package.'+
                            ' Chen2013 and Srinivas2017 are provided with Piperine. [default: Srinivas2017]',
                        type=str)

    parser.add_argument("-x", '--extrapars',
                        help='Parameters sent to SpuriousSSM. Common parameters include "bored", "imax", "tmax" and \
                              can be set by supplying equality statements e.g. "bored=500". [default: None]',
                        type=str)

    parser.add_argument("-q", '--quick',
                        action='store_true',
                        help='Make random numbers instead of computing heuristics to save time[default: False]')
    return parser


def main():
    '''
    Function called by the command line function 'piperine-design'
    '''
    parser = get_design_parser()
    args = parser.parse_args()

    ############## Interpret arguments.
    # Precedence is : command line arguments > crn file arguments > default
    try:
        assert(os.path.isfile(args.crnfile))
    except:
        print('Input file {} does not exist'.format(args.crnfile))
        exit()

    if args.crnfile[-4:] == '.crn':
        basename = args.crnfile[:-4]
        crnfile = args.crnfile
    else:
        basename = args.crnfile.split('.')[0]
        crnfile = args.crnfile[:]

    # Find absolute path to basename
    basedir = os.path.dirname(basename)
    if basedir == '':
        basename = os.getcwd() + os.path.sep + basename
        crnfile = os.getcwd() + os.path.sep + args.crnfile

    # Read parameters from CRN file, only to see if translation scheme is defined
    parameters = get_parameters_form_crn_file(crnfile, None)

    # Import translation scheme package
    if args.translation_scheme:
        translation_scheme = importlib.import_module("."+args.translation_scheme, 'piperine')
    elif 'translation_scheme' in parameters:
        translation_scheme = importlib.import_module("."+parameters['translation_scheme'], 'piperine')
    else:
        translation_scheme = default_translation_scheme

    # metrarank is default
    # Make an error response for providing a method that is not in the approved list
    optimization_methods = ['metarank', 'ranksum']
    if args.optimizer in optimization_methods:
        optimizer = args.optimizer
    else:
        optimizer = 'metarank'

    # Sorry for this horrible line. "translation_scheme" is a package that holds the
    # translation and energetics modules. "translation" is a module that provides the code
    # that generates PIL descriptions of the DNA implementation.
    translation = translation_scheme.translation

    # Make dictionary for design parameters. Start with default parameters
    design_param_dict = dict(zip(translation.param_terms, translation.default_params))

    # Get compilation parameters and define energetics instance.
    # This function uses the translation module to look for translation-specific parameters
    parameters = get_parameters_form_crn_file(crnfile, translation)

    # Apply parameter option preference
    if args.length:
        toehold_length = args.length
    elif 'toehold_length' in parameters:
        toehold_length = parameters['toehold_length']
    elif translation.toehold_length_term in parameters:
        toehold_length = parameters[translation.toehold_length_term]
    else:
        toehold_length = design_param_dict[translation.toehold_length_term]

    design_param_dict[translation.toehold_length_term] = toehold_length

    if args.energy:
        targetdG = args.energy
    elif 'toehold_energy' in parameters:
        targetdG = parameters['toehold_energy']
    else:
        targetdG = 7.7

    if args.deviation:
        deviation = args.deviation
    elif 'toehold_deviation' in parameters:
        deviation = parameters['toehold_deviation']
    else:
        deviation = 0.5

    if args.maxspurious:
        max_spurious = args.maxspurious
    elif 'toehold_spurious' in parameters:
        max_spurious = parameters['toehold_spurious']
    else:
        max_spurious = 0.4

    if args.candidates:
        n = args.candidates
    elif 'n' in parameters:
        n = parameters['n']
    else:
        n = 4

    if args.designparams:
        design_params = args.designparams
    else:
        design_params = translation.default_params
        for term in translation.param_terms:
            if term in parameters:
                design_param_dict[term] = parameters[term]

    try:
        assert(toehold_length == design_param_dict[translation.toehold_length_term])
    except AssertionError:
        raise AssertionError("Toehold length contradictions in input arguments")

    if args.extrapars:
        extra_pars = args.extrapars
    elif 'spurious_design_parameters' in parameters:
        extra_pars = parameters['spurious_design_parameters']
    else:
        extra_pars = ""

    energyfuncs = translation_scheme.energetics.energyfuncs(targetdG=targetdG,
                                         length=toehold_length,
                                         deviation=deviation,
                                         max_spurious=max_spurious)


    out = run_designer(basename,
                     reps=n,
                     design_params=design_params,
                     translation=translation,
                     optimizer=optimizer,
                     energyfuncs=energyfuncs,
                     extra_pars=extra_pars,
                     quick=args.quick)
    winner = out[2]
    if winner is None:
        return None
    else:
        print('Winning sequence set is index {}'.format(winner))
        print('Find this sequence data in the file {}{}.seqs'.format(basename, winner))
        print('The table of scores from this session is in the file {}_scores.csv'.format(basename))
        print('Results of the candidate selection are found in {}_score_report.txt'.format(basename))

def get_scorefixed_parser():
    import argparse
    descr = "Command line utility for scoring a designed sequence set."
    usage = "\n\n\n\
Call template with short option flags. Options are shown in brackets. Capitalized terms stand in for \
required argument or multiple arguments.\n\
piperine-score CRNFILE FIXEDFILE [-e ENERGY] [-p DESIGNPARAMS ...] [-t TRANSLATION_SHEME] [-x EXTRAPARS] [-q] \n\n\
Default execution parameters are stated in the option descriptions. The following are example executions with CRN \
file my_very_own.crn, fixed file my.fixed, and option arguments to override the default settings. \n\
    \n\
    Score the sequence set according to the default translation scheme, Srinivas2017:\n \n\
    piperine-score my_very_own.crn my.fixed\n \
\n \
\n\
    Score another fixed file implemented using the Chen2013 translation scheme:\n\n\
    piperine-score my_very_own.crn chen2013_implementation.fixed -t Chen2013\n\
    or\n\
    piperine-score my_very_own.crn chen2013_implementation.fixed -translation_scheme Chen2013\n\n\n\
    Score the sequence set, but with a different intended toehold binding energy:\n\n\
    piperine-score my_very_own.crn my.fixed -e 9\n\
    "
    parser = argparse.ArgumentParser(description=descr, usage=usage)
    parser.add_argument("crnfile",
                        help='Text file describing CRN. May also define parameters.',
                        type=str)

    parser.add_argument("fixedfile",
                        help='Text file containing sequence constraints. Files read or generated by ' +
                              'this call will have the form: basename.extension ; '+
                          'E.g. basename.fixed may generate basename.pil .',
                        type=str)

    parser.add_argument("-e", "--energy",
                        help='Target toehold binding energy in kcal/mol. [default: 7.7]',
                        type=float)

    parser.add_argument("-p", '--designparams',
                        help='A string of integers that are parameters to the compiling the PIL file from the '+
                             ' system file (.sys extension). [default: Finds in module]',
                        type=int,
                        nargs="*")

    parser.add_argument("-t", '--translation_scheme',
                        help='Provide a string, the name of the Python package describing the translation scheme used to convert '+
                            'the CRN to DNA strands and complexes. See piperine.Srinivas2017 for an example of such a package.'+
                            ' [default: Srinivas2017]',
                        type=str)
    parser.add_argument("-q", '--quick',
                        action='store_true',
                        help='Make random numbers instead of computing heuristics to save time. [default: False]')
    return parser

def score():
    '''
    Function called by the command line function 'piperine-score'
    '''
    parser = get_scorefixed_parser()
    args = parser.parse_args()

    ############## Interpret arguments.
    # Precedence is : command line arguments > crn file arguments > default
    assert(args.crnfile[-4:] == '.crn')
    crnfile = args.crnfile
    fixedfile = args.fixedfile
    basename = fixedfile[:-6]

    # Find absolute path to basename
    basedir = os.path.dirname(basename)
    if basedir == '':
        basename = os.getcwd() + os.path.sep + basename
        crnfile = os.getcwd() + os.path.sep + args.crnfile

    # Read parameters from CRN file, only to see if translation scheme is defined
    parameters = get_parameters_form_crn_file(crnfile, None)

    # Import translation scheme package
    if args.translation_scheme:
        translation_scheme = importlib.import_module("."+args.translation_scheme, 'piperine')
    elif 'translation_scheme' in parameters:
        translation_scheme = importlib.import_module("."+parameters['translation_scheme'], 'piperine')
    else:
        translation_scheme = default_translation_scheme

    # Sorry for this horrible line. "translation_scheme" is a package that holds the
    # translation and energetics modules. "translation" is a module that provides the code
    # that generates PIL descriptions of the DNA implementation.
    translation = translation_scheme.translation

    # Make dictionary for design parameters. Start with default parameters
    design_param_dict = dict(zip(translation.param_terms, translation.default_params))

    # Get compilation parameters and define energetics instance.
    # This function uses the translation module to look for translation-specific parameters
    parameters = get_parameters_form_crn_file(crnfile, translation)

    # Apply parameter option preference
    if args.designparams:
        design_params = args.designparams
    else:
        design_params = translation.default_params
        for term in translation.param_terms:
            if term in parameters:
                design_param_dict[term] = parameters[term]

    if 'toehold_length' in parameters:
        toehold_length = parameters['toehold_length']
    elif translation.toehold_length_term in parameters:
        toehold_length = parameters[translation.toehold_length_term]
    else:
        toehold_length = design_param_dict[translation.toehold_length_term]

    design_param_dict[translation.toehold_length_term] = toehold_length

    if args.energy:
        targetdG = args.energy
    elif 'toehold_energy' in parameters:
        targetdG = parameters['toehold_energy']
    else:
        targetdG = 7.7

    try:
        assert(toehold_length == design_param_dict[translation.toehold_length_term])
    except AssertionError:
        raise AssertionError("Toehold length contradictions in input arguments")

    energyfuncs = translation_scheme.energetics.energyfuncs(targetdG=targetdG,
                                         length=toehold_length)


    out = score_fixed(fixedfile,
                      basename,
                      crn_file=crnfile,
                      design_params=design_params,
                      translation=translation,
                      optimizer=optimizer,
                      energyfuncs=energyfuncs,
                      quick=args.quick)
    print(dict(zip(out[1], out[0])))
