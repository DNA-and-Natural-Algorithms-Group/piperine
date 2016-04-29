# Default parameters
default_params = (7,)
# The .comp files to be imported in the .sys file
comps = ['leakless_and', 'leakless_translate', 'leakless_translate_flux', 'set_leakless_signals_equal']
# A component that sets "signal domains" equal in the form "A->B" for equality
set_equals = 'set_leakless_signals_equal'
# Strings used to generate reaction lines in the .sys file.
# The second item is the argument call
rxn_strings = ['component ', '(<l>)']
# The argument call at the system file header
param_string = '(l): '

# A dictionary mapping number of reactants to information required to specify
# a gate object. 1 leads to nothing, as the Leakless implementation makes
# "react gates" for in the case of AND logic, a 2-input case
react = {1:[],
         2:['leakless_and',                 # comp file
           ['{0}-And', '{0}-And_Waste'],    # Complexes
           ['{0}-F']]                       # Top Strands
        }

# A dictionary mapping the number of reactants to information required
# to specify a produce gate. In the leakless implementation, each reaction
# may only produce one product, and two reactants implies AND logic. There is 
# a consideration to make regarding the input to the translator component. 
# An AND module releases a longer molecule than a signal strand. Even though
# the design of DNA complexes does not need to consider the extra domains of
# the flux molecule, Pepper requires two separate components to achieve
# identical processing of different input strands. 
produ = { 1:
          ['leakless_translate',            # comp file
          ['{0}-T1', '{0}-T1_Waste',        
           '{0}-T2', '{0}-T2_Waste'],       # Complexes
          ['{0}-t1_top']],                   # Top strands
          2:
          ['leakless_translate_flux',       # comp file
          ['{0}-T1', '{0}-T1_Waste',        
           '{0}-T2', '{0}-T2_Waste'],       # Complexes
          ['{0}-t1_top']]}                   # Top strands

# All strands have a pepper_names attribute, which holds a dictionary with 
# the "pepper_keys" as keys and "pepper_values" as values. Generally, a key
# should correspond to something that the Heuristics will reference ( such as
# toeholds, branch migration regions...) or some detail that is rxn scheme
# specific ( eg History domains in the Leakless approach are borrowed from
# reactant strands. The 'history reference' key makes this explicit. The DSD
# approach, however, has unique history domains for each production instance
# of a species, and lacks this key.)
pepper_keys = ['sequence',
               'toeholds',
               'bm domains',
               'history domains',
               'history reference',
               'all domains']

pepper_values = {'r': [['{0}-AB', # Sequence
                        ['{0}-a1', '{0}-b1', '{0}-b3'], # Toeholds
                        ['{0}-b3{0}-b2{0}-b1{0}-a3{0}-a2', 
                         '{0}-b2{0}-b1{0}-a3{0}-a2{0}-a1'], # bm domains
                        [], # hisotry domain
                        ['{0}-b3{0}-b2'], #history reference
                        ['{0}-b3','{0}-b2',
                         '{0}-b1','{0}-a3',
                         '{0}-a2','{0}-a1']], # all domains
                       ['{0}-CD',  # Sequence
                        ['{0}-c1', '{0}-d1', '{0}-d3'], # Toeholds
                        ['{0}-d3{0}-d2{0}-d1{0}-c3{0}-c2', 
                         '{0}-d2{0}-d1{0}-c3{0}-c2{0}-c1'], # bm domains
                        [], # history domain
                        ['{0}-d3{0}-d2'], # history reference
                        ['{0}-d3', '{0}-d2', '{0}-d1',
                         '{0}-c3', '{0}-c2', '{0}-c1']]], # all domains
                 'p': [['{0}-CD',  # Sequence
                        ['{0}-c1', '{0}-d1', '{0}-d3'],# Toeholds
                        ['{0}-d3{0}-d2{0}-d1{0}-c3{0}-c2', 
                         '{0}-d2{0}-d1{0}-c3{0}-c2{0}-c1', 
                         '{0}-c3{0}-c2{0}-c1{0}-b3{0}-b2'],# bm domains
                        ['{0}-b3{0}-b2'],# hisotry domain
                        ['{0}-d3{0}-d2'],
                        ['{0}-d3','{0}-d2','{0}-d1','{0}-c3',
                         '{0}-c2','{0}-c1','{0}-b3','{0}-b2']], 
                       ['{0}-XY', 
                        ['{0}-x1','{0}-y1','{0}-y3'],
                        ['{0}-y3{0}-y2{0}-y1{0}-x3{0}-x2', 
                         '{0}-y2{0}-y1{0}-x3{0}-x2{0}-x1'],
                        ['{0}-c3{0}-c2'],
                        ['{0}-y3{0}-y2'],
                        ['{0}-y3','{0}-y2','{0}-y1','{0}-x3',
                         '{0}-x2','{0}-x1','{0}-c3','{0}-c2']]]}

# Flux strands are a hacky treatment of SignalStrands. Just enough information
# is supplied to allow Flux strands to inherit SignalStrand method attributes.
flux_sequence = '{0}-f'
flux_toehold = '{0}-b1'
# Number of toeholds per signal strand
n_th = 3

def format_list(templates, word):
    if type(templates) is list:
        return [format_list(temp, word) for temp in templates]
    else:
        return templates.format(word)

class SignalStrand(object):
    
    def __init__(self, species):
        self.species = species
        self.history_domains = []
        self.names = []
        self.rxns = []
        self.sequences = []
    
    def __repr__(self):
        base_str = 'species {0} family {1}'
        return base_str.format(self.species, ' '.join(self.names[:]))

    def get_names(self):
        names = self.names[:]
        species = self.species
        names.append(species)
        return names
    
    def set_identity_domains(self, degree, rxn_name, side='r'):
        '''Associates signal identity domains to the object
        
        The reaction processor assigns the first appropriate PIL name to a sequ
        ence. When si gnal strands are encountered, the processor either retrie
        ves PIL names fr om its memory or generates and assigns new PIL names. 
        This function is used in the second case. This function relies on the v
        ariables defined ahead of the Classes.py file to determine the reaction 
        name and .comp file, what the PIL names for each of the signal's domai
        ns will be, and then assigns them to this fledgling species. 

        Args:
            degree: First or second (0 or 1) species in the equation half
            rxn_name: Reaction name appearing in the system file
            side: ('r' or 'p') reactants or products half-equation ('r')
        '''
        self.pepper_names = dict(
            zip(pepper_keys, 
                format_list(pepper_values[side][degree], rxn_name)))
    
    def add_instance(self, name, degree, rxn_name, ref_strand):
        '''Associates a new history domain, sequence, name, and reaction name
        
        When the reaction processor comes upon a signal product, it will
        typically have a unique history domain. This function copies and appends
        information in a reference strand to the history_domains attribute. It
        also saves all available information regarding its location in the sys
        file.
        
        Args:
            name: Name appearing in the system file
            degree: Number of reactants, used to discriminate whether a flux
                    translator or a signal translator is needed
            rxn_name: Reaction name appearing in the system file
            ref_strand: Reactant strand contributing the history domains
        '''
        pepper_names = dict(zip(pepper_keys, 
                                format_list(pepper_values['p'][degree], rxn_name)))
        self.history_domains.append(pepper_names['history reference'])
        self.sequences.append(pepper_names['sequence'])
        self.names.append(name)
        self.rxns.append(rxn_name)
    
    def make_strand_instance(self, name):
        return StrandInstance(self, name)
    
    def th(self, i):
        '''Return references to the requested toehold sequence
        
        Signal strands may use three domains as toeholds. This function
        accepts an index and returns the requested toehold. A signal strand
        will have history domains, then a1 a2 a3 b1 b2 b3 domains, 5->3.
        ['a1', 'b1', 'b3'] is the list drawn from, indexed from 0. 
        Returns:
            * A list of references to toehold sequences
        '''
        return self.pepper_names['toeholds'][i]
    
    def get_bms(self):
        '''Return references to branch migration domains
        
        A branch migration domain is a domain of a top strand that acts to
        attach to a complex "indefinitely", until strand displacement. This
        function returns the references to each branch migration region of the
        strand.
        Returns:
            * A list of references to branch migration domains.
        '''
        return self.pepper_names['bm domains']
    
    def get_top_strands(self):
        '''Return references to each instance of this strand
        
        A signal strand, with associated toeholds and branch migration regions,
        will also have associated history domains from its different produce
        gates. 
        Returns:
            * A list of references to strand instance sequences
        '''
        if self.sequences == []:
            return [self.pepper_names['sequence']]
        return self.sequences[:]
    
    def get_noninteracting_peppernames(self, th, dom_list=None):
        '''
        
        '''
        out_list = []
        for i, hist in enumerate(self.history_domains):
            growing_list = []
            growing_str = hist
            for dom in self.pepper_names['all domains']:
                if th == dom:
                    growing_list.append(growing_str)
                    growing_str = ''
                else:
                    growing_str = growing_str + dom
            out_list.extend(growing_list)
        return out_list

class StrandInstance(SignalStrand):
    def __init__(self, signal, name):
        self.species = signal.species
        self.name = name
        self.names = [name]
        self.pepper_names = signal.pepper_names    
        if name in signal.names:
            i = [ x for x in range(len(signal.names)) if signal.names[x] == name]
            i = i[0]
            self.sequences = signal.sequences[i]
            self.rxns = signal.rxns[i]
            self.history_domains = signal.history_domains[i]

class FluxStrand(SignalStrand):
    def __init__(self, name, rxn_name, ref_strands):
        self.name = name
        self.species = name
        self.names = []
        self.pepper_names = {'sequence':flux_sequence.format(rxn_name)}
        strand_one, strand_two = ref_strands
        ths = [strand_one.th(1), strand_two.th(0)]
        self.pepper_names.update({'toeholds':ths})

class Gate(object):
    '''
    A class to help organize and access pertinent names used by the pepper com
    piler. Children of this class provide all the functionality, this parent
    class is, right now, defined uselessly for future purposes\
    '''
    def __repr__(self):
        return self.get_reaction_line()[:-1]
    
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
        return '{}{} = {}{}: {}\n'.format(rxn_strings[0], rxn, comp, rxn_strings[1], eq)

class ReactGate(Gate):
    def __init__(self, rxn_name, in_strands, out_strands, degree, params):
        self.comp, complexes, top_strands = react[degree]
        self.complexes, self.top_strands = \
            [format_list(x, rxn_name) for x in (complexes, top_strands)]
        self.rxn_name = rxn_name
        self.in_strands = in_strands
        self.out_strands = out_strands
        ############
        l = params[0]
        self.top_s_dict = dict(zip(self.top_strands,[range(1+l*2, 1+l*3+4)]))
        if degree == 2:
            # Grab toeholds, a1 and d3
            self.base = [ in_strands[0].th(0), in_strands[1].th(2) ]
            # Set toe nointeraction map
            self.toe_nointeract_map = dict()
            # First toe of first strand,
            # second toe of second strand do not interact with flux
            self.toe_nointeract_map = dict(
                zip([ in_strands[x].th(x) for x in [0,1]],
                    [self.top_strands[:], self.top_strands[:]]))
        else:
            print 'PROBLEM!!!'

class ProduceGate(Gate):
    def __init__(self, rxn_name, in_strands, out_strands, degree, params):
        self.comp, complexes, top_strands = produ[degree]
        self.complexes, self.top_strands = \
            [format_list(x, rxn_name) for x in (complexes, top_strands)]
        self.rxn_name = rxn_name
        self.in_strands = in_strands
        self.out_strands = out_strands
        ############
        l = params[0]
        self.top_s_dict = dict(zip(self.top_strands,[range(1+l*2, 1+l*3+4)]))
        # Degree refers to number of reactants, dicating whether the
        # translate gate takes in a Flux or a Signal Strand
        if degree == 1:
            strand = in_strands[0]
            # Grab toeholds, b1 and c1
            self.base = [ strand.th(0), strand.th(1) ]
            # When taking in a Signal Strand, the second toehold
            # interacts with the top strand
            self.toe_nointeract_map = dict([ (strand.th(1), [] )])
        elif degree == 2:
            strand = in_strands[0]
            self.base = [ strand.th(0), strand.th(1) ]
            # When taking in a Signal Strand, the second toehold
            # interacts with the top strand
            self.toe_nointeract_map = dict([ (strand.th(1), [] )])

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
    # A dictionary of species to their strand objects
    species_instances = dict()
    # A list to be populated with gate information
    gates = []
    
    # For each reaction, determine the individual react-produce reactions
    # needed to be specified and write these lines to the system file
    for i, rxn in enumerate(rxns):
        I = str(i)
        react_name = 'AND' + I
        produce_name = 'TL' + I
        nr = sum(rxn['stoich_r'])
        if nr == 2:
            reactants = []
            for j, r_species in enumerate(rxn['reactants']):
                if r_species not in species_instances:
                    strand = SignalStrand(r_species)
                    strand.set_identity_domains(j, react_name, 'r')
                    species_instances[r_species] = strand
                else:
                    strand = species_instances[r_species]
                for q in range(rxn['stoich_r'][j]):
                    reactants.append(strand.make_strand_instance(r_species))
                ref = strand
            if nr == 2:
                flux = FluxStrand('Flux'+I, react_name, reactants)
                ref = flux
                gates.append(ReactGate(react_name, reactants, [flux], nr, d_params))
                reactants = [flux]
        else:
            reactants = []
            r_species = rxn['reactants'][0]
            if r_species not in species_instances:
                strand = SignalStrand(r_species)
                strand.set_identity_domains(0, produce_name, 'r')
                species_instances[r_species] = strand
            else:
                strand = species_instances[r_species]
            reactants.append(strand.make_strand_instance(r_species))
            ref = strand
        p_species = rxn['products'][:]
        if rxn['products']:
            p_stoich = rxn['stoich_p']
            if p_stoich == [1]:
                product = p_species[0]
                p_instance = product + I + 'a'
                if product not in species_instances:
                    strand = SignalStrand(product)
                    strand.set_identity_domains(0, produce_name, 'p')
                    species_instances[product] = strand
                else:
                    strand = species_instances[product]
                strand.add_instance(p_instance, nr-1, produce_name, ref)
                products = [strand.make_strand_instance(p_instance)]
                np = 1
            else:
                print 'unexpected stoichiometry!'
                sys.exit('Unexpected product stoichiometry in reaction ' + I)
            gates.append(ProduceGate(produce_name, reactants, products, nr, d_params))
        else:
            # No products! Null output produce gates 
            print 'unexpected stoichiometry!'
            sys.exit('Unexpected product stoichiometry in reaction ' + I)
    
    strands = species_instances.values()
    return (gates, strands)
