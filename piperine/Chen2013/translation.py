import pkg_resources

# comps = [pkg_resources.resource_stream('piperine', 'data/cardelli1-1rxn.comp',
#   'data/cardelli1-2rxn.comp', 'data/cardelli2-1rxn.comp', 'data/cardelli2-2rxn.comp')]
# Default params
default_params = (5, 15)
param_terms = ('t', 'bm')
toehold_length_term = 't'

n_th = 1
# The .comp files to be imported in the .sys file
comps = ['cardelli_1r1p_rxn', 'cardelli_1r2p_rxn',
         'cardelli_2r1p_rxn', 'cardelli_2r2p_rxn']
# Strings used to generate reaction lines in the .sys file.
# The second item is the argument call
rxn_strings = ['component', '(<t>, <bm>)']
# The argument call at the system file header
param_string = '(t, bm): '

pepper_keys = ['sequences',
               'toeholds',
               'bm domains']

dom = ['ta-am', 'tb-bm', 'tc-cm', 'td-dm']
seq = ['a-aR', 'b-bR', 'c-cR', 'd-dR']

def domains(r, p):
    dom_r = dom[:r]
    dom_p = dom[r:r+p][::-1]
    join = "-".join(dom_r) + "-tr-rm-tq"
    fork = "ti-im-" + "-".join(dom_p) + "-tr-rm-tq"
    j = join.split('-')
    f = fork.split('-')
    return (['-'.join(j[i:i+2]) for i in range(len(j)-1)],
            ['-'.join(f[i:i+2]) for i in range(len(f)-2)])

def sequences(r, p):
    seq_r = seq[:r]
    seq_p = seq[r:r+p][::-1]
    return (('-'.join(seq_r)).split('-') + ['helper', 'translator'],
        ['waste', 'end'] + ('-'.join(seq_p)).split('-') + ['helper'])

def format_list(templates, word):
    if type(templates) is list:
        return [format_list(temp, word) for temp in templates]
    else:
        return templates.format(word)

def flatten(in_list):
    if type(in_list) is list or type(in_list) is tuple:
        out_list = []
        for x in in_list:
            out_list.extend(flatten(x))
        return out_list
    else:
        return [in_list]

class SignalStrand(object):

    def __init__(self, species):
        self.species = species
        self.name = species
        self.rxns = []
        self.pepper_names = dict(list(zip(pepper_keys,
                                [[] for _ in range(len(pepper_keys))])))

    def __repr__(self):
        base_str = 'species {0} family {1}'
        return base_str.format(self.species, ' '.join(flatten(self.get_top_strands())))

    def add_instance(self, sequence, domains, rxn_name):
        sequence = rxn_name + '-' + sequence
        toeholds = [rxn_name + '-' + dom for dom in domains.split('-') if 't' in dom]
        branches = [rxn_name + '-' + dom for dom in domains.split('-') if 't' not in dom]
        if not self.pepper_names['sequences']:
            self.pepper_names['sequences'].append(sequence)
        if not self.pepper_names['toeholds']:
            self.pepper_names['toeholds'].extend(toeholds)
        if not self.pepper_names['bm domains']:
            self.pepper_names['bm domains'].extend(branches)
        self.rxns.append(rxn_name)

    def th(self, i):
        return self.pepper_names['toeholds'][i]

    def get_ths(self):
        bm = self.get_top_strands()[0]
        if bm.endswith('R') or bm.endswith('end'): # Ignore reverse and end strands
            return []
        else:
            return self.pepper_names['toeholds'][:]

    def get_bms(self):
        bm = self.get_top_strands()[0]
        if bm.endswith('R') or bm.endswith('end') or bm.endswith('translator'): # Ignore repeats
            return []
        else:
            return self.pepper_names['bm domains'][:]

    def get_top_strands(self):
        return [self.pepper_names['sequences'][0]]

    def get_noninteracting_peppernames(self, th):
        toeholds = self.pepper_names['toeholds']
        bms = self.pepper_names['bm domains']
        if th in toeholds:
            return [t for t in toeholds if t != th] + bms
        else:
            return self.pepper_names['sequences']

def bm_ranges(nr, np, params):
    t, bm = params
    forward = range(1, t + 4)
    reverse = range(bm - 2, bm + t + 1)
    in_ranges = flatten([[forward, reverse] for _ in range(nr + 1)])
    out_ranges = [forward] + flatten([[reverse, forward] for _ in range(np + 1)])
    return [list(el) for el in in_ranges] + [list(el) for el in out_ranges]

def F(strands, rxn_name):
    return {}

class Cardellirxn(object):
    '''
    A class to help organize and access pertinent names used by the pepper com
    piler. Children of this class provide all the functionality, this parent
    class is, right now, defined uselessly for future purposes
    '''
    def __init__(self, rxn_name, in_strands, out_strands, params):
        in_strands, nr = in_strands
        out_strands, np = out_strands
        self.nr = nr
        self.np = np
        self.comp = 'cardelli_{}r{}p_rxn'.format(nr, np)
        self.complexes = format_list(['{0}-Join'] + ['{0}-Join_' + str(i) for i in range(1,nr+2)] +\
                        ['{0}-Fork'] + ['{0}-Fork_' + str(i) for i in range(1,np+3)] +\
                        ['{0}-Reporter', '{0}-Reporter_Waste'], rxn_name)
        self.top_strands = [rxn_name + '-reporter_top']
        self.rxn_name = rxn_name
        self.in_strands = in_strands
        self.out_strands = out_strands
        in_names = [strand.get_top_strands()[0] for strand in in_strands]
        out_names = [strand.get_top_strands()[0] for strand in out_strands]
        # Grab all toeholds
        self.base = ['ta', 'tb', 'tc', 'td'][:nr+np] + ['tr', 'tq', 'ti']
        self.base = [rxn_name + '-' + el for el in self.base]

        # Hard-coded splitting of top-strand domains for toehold occlusion calculation
        self.toe_nointeract_map = F(in_strands + out_strands, rxn_name)
        self.top_s_dict = dict(list(zip(in_names + out_names, bm_ranges(nr, np, params))))

    def __repr__(self):
        return self.get_reaction_line()

    def get_noninteracting_peppernames(self, toehold):
        # Toehold should be the exact name of a toehold domain in the pil/mfe
        # file, which should match the name assigned to the SignalStrand object
        # passed to the gate upon initialization. The is-a Gate object should
        # have an initialized overloading of this function that sees whether
        # the input toehold matches one of the toeholds in the gate, and if so
        # which top strands it should interact with.
        if toehold in self.toe_nointeract_map:
            return self.toe_nointeract_map[toehold]
        else:
            return self.top_strands[:]

    def get_complexes(self):
        return self.complexes[:]

    def get_top_strands(self):
        return self.top_strands[:]

    def get_top_strand_dict(self):
        return self.top_s_dict.copy()

    def get_base_domains(self):
        return self.base[:]

    def get_reaction_line(self):
        comp = self.comp
        rxn = self.rxn_name
        in_names = [self.in_strands[2*i].name for i in range(self.nr)]
        out_names = [self.out_strands[::-1][2*i+2].name for i in range(self.np)]
        eq = ' + '.join(in_names) + ' -> ' + ' + '.join(out_names)
        return '{} {} = {}{}: {}\n'.format(rxn_strings[0], rxn, comp, rxn_strings[1], eq)

def process_rxns(rxns, species, d_params):
    """ Read CRN info into lists of strands and gates

    This function iterates through each reaction recording the gates and
    strands required to model it using DNA. These strand and gate objects
    hold references to specific strands, complexes, domains, and subsequences
    of the DSD system that help construct the inputs to the Pepper compiler
    and spuriousSSM. In addition, they allow the scoring functions to access
    the specific nucleotide sequences they require.

    Args:
        rxns: Output of a import_crn call, a list of reaction dicts
        species: List of unique species in the CRN
        d_params: Parameters to the system
    Returns: A tuple of the following
        gates: A list of gate objects
        strands: A list of strand objects
    """
    import string
    # A list to be populated with gate information
    gates = []

    # A dictionary of species to their strand objects
    species_dict = {}

    # For each reaction, determine the individual react-produce reactions
    # needed to be specified and write these lines to the system file
    for i, rxn in enumerate(rxns):
        I = str(i)
        rxn_name = 'r' + I
        reverse_name = 'Reverse' + I
        fuel_name = 'RXNFuel' + I
        waste_name = 'RXNWaste' + I

        nr = sum(rxn['stoich_r'])
        np = sum(rxn['stoich_p'])
        if nr not in [0, 1, 2] :
            print("Unexpected stoichiometry!")
            raise Exception("Incorrect reactant stochiometry {} for reaction {}. Must be 0, 1, 2".format(nr, I))
        if nr == 0:
            nr += 1
            rxn['stoich_r'].append(1)
            rxn['reactants'].append(fuel_name)
            species.append(fuel_name)
        if np not in [0, 1, 2] :
            print("Unexpected stoichiometry!")
            raise Exception("Incorrect product stochiometry {} for reaction {}. Must be 0, 1, 2".format(np, I))
        if np == 0:
            np += 1
            rxn['stoich_p'].append(1)
            rxn['products'].append(waste_name)
            species.append(waste_name)

        r_species = flatten([rxn['stoich_r'][j]*[spec] for j, spec in enumerate(rxn['reactants'])])
        p_species = flatten([rxn['stoich_p'][j]*[spec] for j, spec in enumerate(rxn['products'])])
        p_species.reverse()
        r_species = flatten(list(zip(r_species, [reverse_name + '_r' + str(j+1) for j in range(len(r_species))])))
        p_species = flatten(list(zip(p_species, [reverse_name + '_p' + str(j+1) for j in range(len(p_species)-1,-1,-1)])))
        r_species = r_species + ['Helper' + I, 'Translator' + I]
        p_species = ['Waste'+I, 'End' + I] + p_species + ['Helper' + I]

        r_sequences, p_sequences = sequences(nr, np)
        r_domains, p_domains = domains(nr, np)

        reactants = []
        for j, spec in enumerate(r_species):
            if spec not in species_dict:
                strand = SignalStrand(spec)
                species_dict[spec] = strand
            else:
                strand = species_dict[spec]
            strand.add_instance(r_sequences[j], r_domains[j], rxn_name)
            reactants.append(strand)

        products = []
        for j, spec in enumerate(p_species):
            if spec not in species_dict:
                strand = SignalStrand(spec)
                species_dict[spec] = strand
            else:
                strand = species_dict[spec]
            strand.add_instance(p_sequences[j], p_domains[j], rxn_name)
            products.append(strand)

        gates.append(Cardellirxn(rxn_name, (reactants, nr), (products, np), d_params))

    strands = list(species_dict.values())
    return (gates, strands)
