# The .comp files to be imported in the .sys file
comps = ['gate_nul_produ', 'gate_uni_produ', 'gate_bim_produ', 
         'gate_uni_react', 'gate_bim_react', 'set_signals_equal']
# A component that sets "signal domains" equal in the form "A->B" for equality
set_equals = 'set_signals_equal'
# Strings used to generate reaction lines in the .sys file.
# The second item is the argument call
rxn_strings = ['component ', '(<t>, <bm>, <c>)']
# The argument call at the system file header
param_string = '(t, bm, c): '

# A dictionary mapping number of reactants to information required to specify
# a gate object: .comp file, refs to complexes, refs to top strands. Indices
# correspond to number of reactants / products
react = {1:['gate_uni_react',
           ['{0}-Gate', '{0}-Gate_waste'],
           ['{0}-Out']],
         2:['gate_bim_react',
           ['{0}-Gate', '{0}-Gate_int', '{0}-Gate_waste'],
           ['{0}-Out', '{0}-Backward']]}

# A dictionary mapping the number of reactants to information required
# to specify a produce gate. 
produ = [[ 'gate_nul_produ', # comp file
          ['{0}-Trans', '{0}-Trans_waste'],      # Complexes 
          ['{0}-C']],                # Top strands
         [ 'gate_uni_produ',           
          ['{0}-Trans', '{0}-Trans_int',                   
           '{0}-Trans_waste','{0}-Trans_cat_waste'],    
          ['{0}-cat_helper', '{0}-helper', '{0}-d']], 
         [ 'gate_bim_produ',                         
          ['{0}-Trans', '{0}-Trans_int',             
           '{0}-Trans_waste','{0}-Trans_cat_waste'], 
          ['{0}-cat_helper', '{0}-helper']]]

pepper_keys = ['sequence',
               'toeholds',
               'bm domains',
               'history domains']

pepper_values = {'r': [['{0}-a', 
                        ['{0}-toe-fa', '{0}-toe-sa'],
                        ['{0}-am{0}-cam'],
                        []], 
                       ['{0}-b', 
                        ['{0}-toe-fb', '{0}-toe-sb'],
                        ['{0}-bm{0}-cbm'],
                        []]], 
                 'p': [['{0}-C', 
                        ['{0}-toe-fc', '{0}-toe-sc'], 
                        ['{0}-cm{0}-ccm'], 
                        ['{0}-ch{0}-cch']],
                       ['{0}-D', 
                        ['{0}-toe-fd', '{0}-toe-sd'], 
                        ['{0}-dm{0}-cdm'],
                        ['{0}-dh{0}-cdh']]]}

# Flux strands are a hacky treatment of SignalStrands. Just enough information
# is supplied to allow Flux strands to inherit SignalStrand method attributes.
flux_sequence = '{}-Out'
flux_toehold = '{}-toe-sb'
# Number of toeholds per signal strand
n_th = 2

def format_list(templates, word):
    if type(templates) is list:
        return [format_list(temp, word) for temp in templates]
    else:
        return templates.format(word)

class SignalStrand(object):
    
    def __init__(self, species):
        self.species = species
        self.names = []
        self.rxns = []
        self.sequences = []
    
    def __repr__(self):
        base_str = 'species {0} family {1}'
        return base_str.format(self.species, ' '.join(self.names[:]))
    
    def set_identity_domains(self, degree, rxn_name, gate='r'):
        self.pepper_names = dict(zip(pepper_keys, 
                                format_list(pepper_values[gate][degree], 
                                            rxn_name)))
    
    def add_instance(self, name, degree, rxn_name):
        pepper_names = dict(zip(pepper_keys, 
                                format_list(pepper_values['p'][degree], 
                                            rxn_name)))
        hd = pepper_names['history domains']
        self.pepper_names['history domains'].append(hd)
        self.sequences.append(pepper_names['sequence'])
        self.names.append(name)
        self.rxns.append(rxn_name)
    
    def make_strand_instance(self, name):
        return StrandInstance(self, name)
    
    def th(self, i):
        return self.pepper_names['toeholds'][i]
    
    def get_names(self):
        names = self.names[:]
        species = self.species
        names.append(species)
        return names
    
    def get_bms(self):
        return self.pepper_names['bm domains']
    
    def get_top_strands(self):
        if len(self.sequences) == 0:
            return [self.pepper_names['sequence']]
        return self.sequences[:]
    
    def get_noninteracting_peppernames(self, th):
        toeholds = self.pepper_names['toeholds']
        bms = self.pepper_names['bm domains']
        hds = self.pepper_names['history domains']
        outlist = []
        for i, hd in enumerate(hds):
            if th == toeholds[0]:
                outlist.extend([hd, toeholds[1] + bm[0]])
            if th == toeholds[1]:
                outlist.extend([bm[0] + toeholds[0] + hd])
            else:
                outlist.append(self.sequences[i])
        return outlist

class StrandInstance(SignalStrand):
    def __init__(self, signal, name):
        self.species = signal.species
        self.name = name
        self.names = [name]
        self.pepper_names = signal.pepper_names.copy()
        if name in signal.names:
            i = [ x for x in range(len(signal.names)) if signal.names[x] == name]
            i = i[0]
            self.sequences = signal.sequences[i]
            self.rxns = signal.rxns[i]

class FluxStrand(SignalStrand):
    def __init__(self, name, rxn_name, ref_strand):
        self.names = [name]
        self.name = name
        self.species = 'Flux'
        self.pepper_names = {'sequence':flux_sequence.format(rxn_name),
                             'toeholds':[ref_strand.th(1)]}

class Gate(object):
    '''
    A class to help organize and access pertinent names used by the pepper com
    piler. Children of this class provide all the functionality, this parent
    class is, right now, defined uselessly for future purposes\
    '''
    
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
        if len(self.in_strands) == 2:
            eq = '{0} + {1}'.format(self.in_strands[0].name,
                                    self.in_strands[1].name)
        else:
            eq = self.in_strands[0].name
        eq = eq + ' -> '
        if len(self.out_strands) == 2:
            eq = eq + '{0} + {1}'.format(self.out_strands[0].name,
                                         self.out_strands[1].name)
        elif len(self.out_strands) == 1:
            eq = eq + self.out_strands[0].name
        return '{} {} = {}{}: {}\n'.format(rxn_strings[0], rxn, comp, rxn_strings[1], eq)

class ReactGate(Gate):
    def __init__(self, rxn_name, in_strands, out_strands, degree, params):
        self.comp, complexes, top_strands = react[degree]
        self.complexes, self.top_strands = \
            [format_list(x, rxn_name) for x in (complexes, top_strands)]
        self.rxn_name = rxn_name
        self.in_strands = in_strands
        self.out_strands = out_strands
        #############
        t, bm, c = params
        if degree == 1:
            self.top_s_dict = dict(zip(self.top_strands,[range(1+bm, 1+bm+t+4)]))
        else:
            self.top_s_dict = dict(zip(self.top_strands,
                                        [range(1+bm, 1+bm+t+4),
                                         range(1+bm+t-2, 1+bm+2*t)]))
        if degree == 1:
            # Grab toeholds, first toehold of first strand
            self.base = [ in_strands[0].th(0)]
            # second toehold interacts with flux
            strand = in_strands[0]
            self.toe_nointeract_map = dict([ (strand.th(1), [] )])
            # Specify complexes
        elif degree == 2:
            # Grab toeholds, first toehold of both strands
            self.base = [ in_strands[0].th(0), in_strands[1].th(0) ]
            # Set toe nointeraction map
            self.toe_nointeract_map = dict()
            # toe-fa does not interact with either top strand
            # toe-sa does not interact with flux
            # toe-fb does not interact with flux
            # toe-sb does not interact with backwards
            self.toe_nointeract_map = dict(
                zip([ in_strands[x].th(y) for x in [0,1] for y in [0, 1]],
                    [self.top_strands[:], self.top_strands[1:], 
                     self.top_strands[1:], self.top_strands[:1]]))

class ProduceGate(Gate):
    def __init__(self, rxn_name, in_strands, out_strands, degree, params):
        self.comp, complexes, top_strands = produ[degree]
        self.complexes, self.top_strands = \
            [format_list(x, rxn_name) for x in (complexes, top_strands)]
        self.rxn_name = rxn_name
        self.in_strands = in_strands
        self.out_strands = out_strands
        #############
        t, bm, c = params
        if degree == 0:
            self.top_s_dict = dict(zip(self.top_strands, [range(1+bm, 1+bm+t+4)]))
        elif degree == 1:
            self.top_s_dict = dict(zip(self.top_strands, [range(1+bm, 1+bm+t+4),
                                                range(1+bm, 1+bm+t+4),
                                                []]))
        else:
            self.top_s_dict = dict(zip(self.top_strands, [range(1+bm, 1+bm+t+4),
                                                range(1+bm, 1+bm+t+4)]))
        self.base = []
        if degree == 0:
            # This strand has no toehold, so should not interact with anything.
            # Therefore, no toehold should appear in the dicionary keys, meaning
            # there are no exceptions to the no-interaction assumption
            self.toe_nointeract_map = dict()
        elif degree == 1:
            self.base = [ in_strands[0].th(0)]
            # With one output, the two toeholds are the flux toehold and first
            # output toehold. The flux, however, is counted elsewhere. So,
            # the first output toehold should interact with all top strands
            strand = out_strands[0]
            self.toe_nointeract_map = dict([ (strand.th(0), [] )])
        elif degree == 2:
            self.base = [ in_strands[0].th(0), out_strands[0].th(0)]
            # Both first toeholds interact with both top strands
            self.toe_nointeract_map = dict(
                        zip( [ x.th(0) for x in out_strands],
                             [ [], [] ] ) )

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
    # A dictionary of species to their strand objects
    species_instances = dict()
    # A list to be populated with gate information
    gates = []
    
    # For each reaction, determine the individual react-produce reactions
    # needed to be specified and write these lines to the system file
    for i, rxn in enumerate(rxns):
        I = str(i)
        react_name = 'r' + I
        produce_name = 'p' + I
        if rxn['stoich_r'] == []:
            nr = 1
            r_species = 'Auxilliary'
            p_species = 'Auxilliary' + I + 'b'
            aux_strand = SignalStrand(p_species)
            aux_strand.set_identity_domains(0, react_name)
            species_instances.update({p_species:aux_strand})
            reactants = [aux_strand.make_strand_instance(p_species)]
            ref_strand = aux_strand
            rxn['stoich_r'] = [1]
            rxn['stoich_p'].append(1)
            rxn['reactants'].append(p_species)
            rxn['products'].append(r_species)
            ref_strand = aux_strand
        elif len(rxn['stoich_r']) <= 2:
            nr = sum(rxn['stoich_r'])
            reactants = []
            for j, r_species in enumerate(rxn['reactants']):
                if r_species not in species_instances:
                    strand = SignalStrand(r_species)
                    strand.set_identity_domains(2 - nr + j, react_name)
                    species_instances[r_species] = strand
                else:
                    strand = species_instances[r_species]
                for q in range(rxn['stoich_r'][j]):
                    reactants.append(strand.make_strand_instance(r_species))
            ref_strand = reactants[nr-1]
        else:
            print 'unexpected stoichiometry!'
            sys.exit('Unexpected reactant stoichiometry in reaction ' + I)
        flux = FluxStrand('Flux'+I, react_name, ref_strand)
        gates.append(ReactGate(react_name, reactants, [flux], nr, d_params))
        
        if rxn['products']:
            suffixes = [I+l for l in string.ascii_lowercase]
            stoichs = rxn['stoich_p']
            prods = rxn['products']
            if stoichs not in [ [1], [1,1], [2] ]:
                stoichs = [ 1 for x in prods ]
            products = []
            np = 0
            for j, stoich, product in zip(range(len(prods)), 
                                                  stoichs, 
                                                  prods):
                if product not in species_instances:
                    strand = SignalStrand(product)
                    strand.set_identity_domains(j, produce_name, 'p')
                    species_instances[product] = strand
                else:
                    strand = species_instances[product]
                for q in range(stoich):
                    instance = product + suffixes[j]
                    strand.add_instance(instance, j, produce_name)
                    products.append(strand.make_strand_instance(instance))
                    np = np + 1
        else:
            # No products! Null output produce gates 
            products = []
            np = 0
        gates.append(ProduceGate(produce_name, [flux], products, np, d_params))
    strands = species_instances.values()
    return (gates, strands)
