from __future__ import division

import numpy as np
import re
from time import time
import stickydesign as sd
import energyfuncs_james as efj

def flatten(x):
    if type(x) in [list, tuple]:
        return [ z for y in x for z in flatten(y)]
    else:
        return [x]
    
    
def get_toeholds(n_ths=6, thold_l=int(7.7), thold_e=7.7, e_dev=0.5, m_spurious=0.4, 
                 e_module=efj, timeout=8):
    """ Generate specified stickyends for the Soloveichik DSD approach
    
    A given run of stickydesign may not generate toeholds that match the 
    Soloveichik approach. This function reports whether the run was successful
    or not. Otherwise, the toeholds are matched to respect the backwards 
    strand back-to-back toeholds. 
    
    Args:
        ef: Stickydesign easyends object
        n_ths: Number of signal strand species to generate toeholds for
        thold_l: Nucleotides in the toeholds
        thold_e: Target binding energy for the toeholds in kcal/mol
        e_dev: Allowable energy deviation in kcal/mol
        m_spurious: Maximum spurious interaction strength as a ratio of target
                    energy
    Returns:
        ends_all: Stickydesign stickyends object
        (e_avg, e_rng): Average and range (max minus min) of toehold energies
    """
    # Give the energetics instance the target energy
    ef = e_module.energyfuncs(targetdG = thold_e)

    # Give StickyDesign a set of trivial, single-nucleotide toeholds to avoid poor
    # designs. I'm not sure if this helps now, but it did once.
    avoid_list = [i * int(thold_l + 2) for i in ['a', 'c', 't']]
    
    # Generate toeholds
    notoes = True
    startime = time()
    while notoes:
        try:
            ends = sd.easyends('TD', thold_l, interaction=thold_e,
                               fdev=e_dev/thold_e, alphabet='h', adjs=['c','g'],
                               maxspurious=m_spurious, energetics=ef,
                               oldends=avoid_list)
            notoes = len(ends) < n_ths
        except ValueError as e:
            if (time() - startime) > timeout:
                return -1
            noetoes = True
    
    th_cands = ends.tolist()
    # remove "avoid" sequences
    th_cands = th_cands[len(avoid_list):]
    # Make as many end in c as possible
    th_cands = th_cands[:n_ths]
    th_all = [ th[1:-1] for th in th_cands]
    ends_full = sd.endarray(th_cands, 'TD')
    ends_all = sd.endarray(th_all, 'TD')
    return ends_all.tolist()

def score_toeholds(toeholds, targetdG=7.7, e_module=efj):
    toeholds_flanked = [ 'c' + th.lower() + 'c' for th in toeholds]
    ef = e_module.energyfuncs(targetdG=targetdG)
    ends = sd.endarray(toeholds_flanked, 'TD')
    e_vec = ef.matching_uniform(ends)
    e_vec_ext = ef.th_external_dG(ends)
    e_vec_int = ef.th_internal_dG(ends)
    e_vec_avg = (e_vec_ext + e_vec_int) / 2
    e_vec_all = np.concatenate( (e_vec_int, e_vec_ext))
    e_avg = e_vec_all.mean()
    e_rng = e_vec_all.max() - e_vec_all.min()
    return (e_avg, e_rng)
    

def generate_pairs(n_ths=3, thold_l=int(7), thold_e=7.7, e_dev=0.5, m_spurious=0.4, 
                 e_module='energyfuncs_james', labels=None):
    """ Generate specified stickyends for the Soloveichik DSD approach
    
    A given run of stickydesign may not generate toeholds that match the 
    Soloveichik approach. This function calls new_toeholds until appropriate
    toeholds are produced
    
    Args:
        n_ths: Number of signal strand species to generate toeholds for
        thold_l: Nucleotides in the toeholds
        thold_e: Target binding energy for the toeholds in kcal/mol
        e_dev: Allowable energy deviation in kcal/mol
        m_spurious: Maximum spurious interaction strength as a ratio of target
                    energy
        e_module: String name of the module to import for energetics functions
        labels: Species names to be passed to the barplot maker
    Returns:
        zip(first_ths, second_ths): List of pair-tuples of matched toeholds
        scores: Average and range (max minus min) of toehold energies
    """
    exec('import ' + e_module + ' as energetics')
    # Give the energetics instance the target energy
    ef = energetics.energyfuncs()
    ef.targetdG = thold_e
    
    # new_toeholds returns 1 when second toehold assignment fails. Use a while 
    # loop to run the toehold generation until successful. 
    ends_all = None
    while type(ends_all) is not sd.endarray:
        output = new_toeholds(ef, n_ths, thold_l, thold_e, e_dev, \
                              m_spurious, labels=labels)
        if type(output) is tuple:
            ends_all, scores = output
        else:
            if output == 1:
                print 'Not enough appropriate first toeholds found. Trying again.'
            if output == 2:
                print 'First-second toehold pairing failed. Trying again.'
            if output == 3:
                print "Dunno"
    ends_list = ends_all.tolist()
    first_ths = ends_list[0::2]
    second_ths = ends_list[1::2]
    return zip(first_ths, second_ths), scores

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-l", "--length", help='Toehold length.[7]', type=int)
    parser.add_argument("-e", "--energy", help='Target toehold binding energy'+
                        'in kcal/mol.[7.7]', type=float)
    parser.add_argument("-d", "--deviation", help='Toehold binding energy st'+\
                        'andard deviation limit in kcal/mol.[0.5]', type=float)
    parser.add_argument("-s", "--maxspurious", help='Maximum spurious interac'+
                        'tion energy as a multiple of target binding energy.['+
                        '0.4]', type=float)
    parser.add_argument("-g", '--energetics', help='The name of the energetics'+
                        ' module to be used for toehold generation.[energyfunc'+
                        's_james]', type=str)
    parser.add_argument("-b", '--barplot', help='User-specified name for' +\
                        ' barplot image file [binding_energy.png]', type=str)
    parser.add_argument("--h", '--heatmap', help='User-specified name for' +\
                        ' spurious heatmap image file [spurious.png]', type=str)
    parser.add_argument("-n", "--number", help='Toeholds needed.[3]', type=int)
    args = parser.parse_args()
    
    ############## Interpret arguments
    # Set toehold length in base pairs
    if args.length:
        thold_l = int(args.length)
        print 'Using length {0}'.format(thold_l)
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
        eval('import ' + args.energetics + 'as ef')
    else:
        import energyfuncs_james as energetics
    
    # Set user-specified barplot image file name
    if args.barplot:
        bplotfile = args.barplot
    else:
        bplotfile = 'binding_energy.png'
    
    # Set user-specified spurious interaction heatmap image file name
    if args.barplot:
        hmfile = args.heatmap
    else:
        hmfile = 'spurious.png'
    
    # Set user-specified spurious interaction heatmap image file name
    if args.number:
        n_ths = args.number
    else:
        n_ths = 3
    
    ############## end Interpret arguments
    
    # Give the energetics instance the target energy
    ef = energetics.energyfuncs()
    ef.targetdG = thold_e
    
    # new_toeholds returns 1 when second toehold assignment fails. Use a while 
    # loop to run the toehold generation until successful. 
    ends_all = None
    while type(ends_all) is not sd.endarray:
        ends_all = new_toeholds(ef, n_ths, thold_l, thold_e, e_dev, m_spurious)
        if type(ends_all) == int:
            if ends_all == 1:
                print 'Not enough appropriate first toeholds found. Trying again.'
            if ends_all == 2:
                print 'First-second toehold pairing failed. Trying again.'
            if ends_all == 3:
                print "dunno"
        
    # Plot energy array
    e_array = sd.energy_array_uniform(ends_all, ef)
    
    print 'Done!'
