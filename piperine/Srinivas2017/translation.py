from __future__ import division, print_function
import pkg_resources

# Default design parameters
default_params = (7, 15, 2)
param_terms = ('t', 'bm', 'c')
toehold_length_term = 't'

# The .comp files to be imported in the .sys file
comps = ['bimrxn']

# Strings used to generate reaction lines in the .sys file.

# The second item is the argument call
param_str = '(<t>, <bm>, <c>)'

# The argument call at the system file header
param_string = '(t, bm, c): '

# List containing component name, list of complexes, and list of strands
pepper_templates =  ['bimrxn',
                    ['{0}-Gate', '{0}-Gate_int', '{0}-Gate_waste',
                     '{0}-Trans', '{0}-Trans_int', '{0}-Trans_waste','{0}-Trans_cat_waste'],
                    ['{0}-Out', '{0}-Backward', '{0}-cat_helper', '{0}-helper']]

# Gate and Strand objects store PIL names of important domains. These domains
# are held in a dictionary with the following keys and values.
pepper_keys = ['sequence',
               'toeholds',
               'bm domains',
               'history domains']

# String templates for important PIL domain names. This is a list of list,
# where the items of the outermost list correspond, in order, to the following
# strands of the bimolecular reaction: first input, second input, first output,
# second output. Each of these members are a list that hold templates for PIL
# names of the following: full strand, list of toehold domains, branch
# migration domains, and history domains.
pepper_values = [['{0}-a',
                  ['{0}-toe-fa', '{0}-toe-sa'],
                  ['{0}-am{0}-cam'],
                  []],
                 ['{0}-b',
                  ['{0}-toe-fb', '{0}-toe-sb'],
                  ['{0}-bm{0}-cbm'],
                  []],
                 ['{0}-c',
                  ['{0}-toe-fc', '{0}-toe-sc'],
                  ['{0}-cm{0}-ccm'],
                  ['{0}-ch{0}-cch']],
                 ['{0}-d',
                  ['{0}-toe-fd', '{0}-toe-sd'],
                  ['{0}-dm{0}-cdm'],
                  ['{0}-dh{0}-cdh']]]


# Inputs to the Toehold Occlusion heuristic function include sequences that
# do not correspond to domains in the component file. These subsequences are
# generated from the string templates stored in the variable specific_names.
# The strings in specific_names are concatenations of PIL domain names that
# match strand definitions in the component file. The input to the heuristic
# score is produced by in-place-replacement of each PIL domain name with its
# nucleotide sequence.  The function F, defined below, produces lists of
# concatenated PIL names. Each list is associated with a toehold domain such
# that the concatenated PIL names in a list do not include the list's
# associated toehold domain. F does this by splitting the concatenated domain
# names into two strings, dividing where the associated toehold domain name
# occurs. The toehold domain name is thus removed.
# The variable nicknames is a list of PIL names for strands, just as the
# variable specific_names is a list of concatenated PIL domain names. In fact,
# both variables define the same strands in the same order. Whenever F finds
# that it does not need to split a concatenated PIL name, it replaces that
# string with the string from the variable nicknames that describes the same
# strand.
nicknames = ['{0}-Out', '{0}-Backward', '{0}-cat_helper', '{0}-helper']
specific_names = ['{0}-ch{0}-cch{0}-toe-sb{0}-bm{0}-cbm',
                  '{0}-toe-fb-suffix{0}-toe-sa{0}-am{0}-cam',
                  '{0}-toe-fd{0}-dh{0}-cdh{0}-toe-fc{0}-ch{0}-cch',
                  '{0}-toe-fd{0}-dh{0}-cdh{0}-toe-fc']
toeholds = [['{0}-toe-fa'], ['{0}-toe-sa'], ['{0}-toe-fb-suffix', '{0}-toe-fb'], ['{0}-toe-sb'],
            ['{0}-toe-fc'], ['{0}-toe-sc'], ['{0}-toe-fd'], ['{0}-toe-sd']]

def format_list(templates, word):
    ''' Recursive function to perform a simple string formatting for lists of lists of strings

        For nested lists containing only strings, this function maintains the
        nesting structure of the input templates doing a single-item format of the
        member strings using argument word.

        Arguments
            templates : Nested list of strings that all have {0} terms for replacement
            word : Substitution string
        Returns
            nested lists of formatted strings
    '''
    if type(templates) is list:
        return [format_list(temp, word) for temp in templates]
    else:
        return templates.format(word)

def flatten(in_list):
    ''' Turns a nested list into a single, un-nested list
    '''
    if type(in_list) is list:
        out_list = []
        for x in in_list:
            out_list.extend(flatten(x))
        return out_list
    else:
        return [in_list]

def F(ordered_species, rxn_name):
    ''' Prepares input to Toehold Occulsion heuristic function.

        Toehold Occlusion calculation requires a dictionary associating a
        toehold with a list of subsequences of top strands that do not contain
        that toehold. However, this top strand list cannot be generated through
        nucleotide sequence comparison because unintentional and intentional
        toehold sequences in a top strand would be indistinguishable. Rather, F
        does this by comparing PIL domain names.

        This function uses variables defined outside of the function body to
        construct the string templates for toehold and strand definitions.

        Arguments
            ordered_species: Strand objects corresponding to signal species
                             associated with a bimolecular reaction.  in the
                             following order:
                                0: First input
                                1: Second input
                                2: First output
                                3: Second output
            rxn_name: String that is entered into the {0} tokens in the strings.
        Returns
            A dictionary with toehold PIL names as keys and a list of
            concatenated PIL domains that correspond to top strand subesquences
            that do not contain the associated toehold domain.
    '''
    # Format the PIL name templates with the reaction name
    th_names = format_list(toeholds, rxn_name)
    domains = format_list(specific_names, rxn_name)
    strands = format_list(nicknames, rxn_name)

    # This loop associates the PIL names for toehold domains from SignalStrand
    # objects to the PIL names for the identical toehold domains that appear in
    # the reaction specified by the argument rxn_name.
    toehold_splits_map = {}
    for i, spec in enumerate(ordered_species):
        spec_first_th = spec.th(0)
        # make new dictionary entry for new strand
        if spec_first_th not in toehold_splits_map:
            toehold_splits_map.update({spec.th(0):th_names[2*i],
                                       spec.th(1):th_names[2*i+1]})
        # update existing dictionary entry
        else:
            toehold_splits_map[spec.th(0)].extend(th_names[2*i])
            toehold_splits_map[spec.th(1)].extend(th_names[2*i+1])

    # The recursive function and the helper function perform the string splitting
    def split_recurse(split_list, re_domains, i):
        ''' Recursive function that splits terms in re_domains at terms found in split_list.

            Arguments
                split_list: A list of PIL names as strings.
                re_domains: A nested list of concatenated PIL names
                i: Top strand index. When no splits are necessary, this index
                   is used to retrieve the strand nickname.
            Returns
                List of concatenated PIL names that do not contain the PIL
                names in split_list.
        '''
        if len(split_list) == 0:
            return re_domains
        if type(split_list) is list:
            split_term = split_list[0]
            splitted = [split_recurse(split_list[1:], x.split(split_term), i) for x in re_domains]
            splitted = flatten(splitted)
            splitted = [x for x in splitted if len(x)>0]
            if len(splitted)>0 and splitted[0] in domains:
                return strands[i]
            return splitted

    def split_helper(split_list):
        ''' Start recursive splitting.
        '''
        results = [ split_recurse(split_list[:], [domains[i]], i) for i in range(len(domains)) ]
        return flatten(results)

    toehold_no_interact_map = {}
    for key in toehold_splits_map:
        toehold_no_interact_map.update({key:split_helper(toehold_splits_map[key])})

    return toehold_no_interact_map

class SignalStrand(object):
    ''' One instance of this class is constructed for each unique signal species
        in the input CRN. The data attributes are assigend values before
        sequence generation and the method attributes procure this data when
        accumulating inputs to the heuristic scoring functions.

        The SignalStrand class was created to
    '''
    def __init__(self, species):
        self.species = species
        self.name = species
        self.rxns = []
        self.sequences = []
        self.pepper_names = {}

    def __repr__(self):
        base_str = 'species {0} family {1}'
        return base_str.format(self.species, ' '.join(self.sequences[:]))

    def set_identity_domains(self, degree, rxn_name):
        ''' When a strand is first encountered, assign PIL names to class
            attributes representing branch migration and toehold domains.

            Arguments
                self: Class instance (automatically supplied).
                degree: Index representing the signal species' position in the
                        bimolecular reaction.
                rxn_name: Reaction prefix as seen in the .sys file.
            Returns
                None
        '''
        self.pepper_names = dict(list(zip(pepper_keys,
                                format_list(pepper_values[degree],
                                            rxn_name))))

    def add_instance(self, degree, rxn_name):
        ''' Records history domains, which are unique to each reaction product.

            Arguments
                self: Class instance (automatically supplied).
                degree: Index representing the signal species' position in the
                        bimolecular reaction.
                rxn_name: Reaction prefix as seen in the .sys file.
            Returns
                None
        '''
        pepper_names = dict(list(zip(pepper_keys,
                                format_list(pepper_values[degree],
                                            rxn_name))))
        hd = flatten(pepper_names['history domains'])
        if hd[0] not in self.pepper_names['history domains']:
            self.pepper_names['history domains'].extend(hd)
        self.sequences.append(pepper_names['sequence'])
        self.rxns.append(rxn_name)

    def th(self, i):
        ''' Returns PIL name for toehold domains.

            Arguments
                self: Class instance (automatically supplied).
                i: Toehold index queried.
            Returns
                PIL name for requested toehold domain, as string.
        '''
        return self.pepper_names['toeholds'][i]

    def get_ths(self):
        ''' Returns list of all PIL names for toehold domains.

            Arguments
                self: Class instance (automatically supplied).
            Returns
                List of all PIL names for toehold domains.
        '''
        return self.pepper_names['toeholds'][:]

    def get_bms(self):
        ''' Returns list of all PIL names for branch migration domains.

            Arguments
                self: Class instance (automatically supplied).
            Returns
                List of all PIL names for branch migration domains.
        '''
        return self.pepper_names['bm domains'] + self.pepper_names['history domains']

    def get_top_strands(self):
        ''' Returns list of all PIL names for all unique strands representing
            this abstract species.

            Arguments
                self: Class instance (automatically supplied).
            Returns
                List of all PIL names for all unique strands representing
                this abstract species.
        '''
        if len(self.sequences) == 0:
            return [self.pepper_names['sequence']]
        return self.sequences[:]

    def get_noninteracting_peppernames(self, th):
        ''' Returns a list of concatenated PIL names as strings that do not
            contain the passed string.

            Rather than employ F, which is used by the get_noninteracting_peppernames
            method of the Bimrxn class and relies on global variables, this function
            assembles concatenated PIL names from the data stored in the pepper_names
            class attribute. When finding the subsequences of top strands that do not
            include a given toehold domain, there are three options:

                1) The signal strand does not have the toehold domain.
                2) The given toehold domain is the first toehold of the signal strand.
                3) The given toehold domain is the second toehold of the signal strand.

            This function returns a list of concatenated PIL names, one for each
            history domain associated with a signal species.


            Arguments
                self: Class instance (automatically supplied).
                th: String. Intended be PIL names for toehold domains taken from
                    either .get_th or .th class methods of SignalStrand instances.
            Returns
                A list of strings, the concatenated PIL names for top strand
                subsequences that do not contain the input toehold domain.
        '''
        # Retrieve relevent PIL names
        toeholds = self.pepper_names['toeholds']
        bms = self.pepper_names['bm domains']
        hds = flatten(self.pepper_names['history domains'])
        outlist = []
        if len(hds) == 0:
            if th == toeholds[0]:
                outlist.append(toeholds[1] + bms[0])
            elif th == toeholds[1]:
                outlist.append(bms[0] + toeholds[0])
            else:
                outlist.append(self.pepper_names['sequence'])
            return outlist
        for i, hd in enumerate(hds):
            if th == toeholds[0]:
                outlist.extend([hd, toeholds[1] + bms[0]])
            elif th == toeholds[1]:
                outlist.append(bms[0] + toeholds[0] + hd)
            else:
                if len(self.sequences) == 0:
                    outlist.append(self.pepper_names['sequence'])
                else:
                    outlist.append(self.sequences[i])
        return outlist

class Bimrxn(object):
    '''
    A class to help organize and access pertinent names used by the pepper com
    piler. Children of this class provide all the functionality, this parent
    class is, right now, defined uselessly for future purposes
    '''
    def __init__(self, rxn_name, in_strands, out_strands, params):
        self.comp, complexes, top_strands = pepper_templates
        self.complexes, self.top_strands = \
            [format_list(x, rxn_name) for x in (complexes, top_strands)]
        self.rxn_name = rxn_name
        self.in_strands = in_strands
        self.out_strands = out_strands
        #############
        t, bm, c = params
        self.params = params
        # Grab all first toeholds, second toehold of second input
        self.base = flatten([[ in_strands[x].th(y), out_strands[x].th(y)] for x in [0,1] for y in [0]] +
                    [ in_strands[1].th(1)])
        # self.base = flatten([ out_strands[x].th(y) for x in [0,1] for y in [0]] +
        #                     [ in_strands[0].th(0), in_strands[1].th(0) + "-suffix"] +
        #                       [ in_strands[1].th(1)])
        # Hard-coded splitting of top-strand domains for toehold occlusion calculation
        self.toe_nointeract_map = F(in_strands + out_strands, rxn_name)
        self.top_s_dict = dict(list(zip(self.top_strands,
                                    [list(range(bm-2, bm+t+1)),
                                     list(range(1, t+3)), # Shorter due to truncated toehold
                                     list(range(t+bm-2, t*2+bm+4)),
                                     list(range(t+bm-2, t*2+bm+1))])))

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

    def get_naive_complexes(self):
        react_gate = 'name+{0}-gate_base+{0}-out+{0}-backward'.format(self.rxn_name)
        react_names = 'React_{1}{2}{0}-c+{0}-react_base+flux_{2}{0}-c+back_{1}{2}'
        react_names = react_names.format(self.rxn_name,self.in_strands[0].name,self.in_strands[1].name)
        produce_gate = 'name+{0}-trans_base+{0}-d+{0}-c+{0}-helper+{0}-cat_helper'.format(self.rxn_name)
        produce_names = 'Produce_{1}{0}-c{0}-d+{0}-produce_base+{0}-d+{0}-c+helper_{2}{0}-d+cat_helper_{0}-c{0}-d'
        produce_names = produce_names.format(self.rxn_name, self.in_strands[1].name, self.out_strands[0].name)
        for strand in self.out_strands:
            strand_versions = strand.get_top_strands()
            for i, ver in enumerate(strand_versions):
                strand_name = "{}{}".format(strand.name, i)
                produce_names = produce_names.replace(ver, strand_name)
                react_names = react_names.replace(ver, strand_name)

        react_terms = [x.split('+') for x in [react_names, react_gate]]
        react_dict = dict(zip(react_terms[0][1:], react_terms[1][1:]))

        produce_terms = [x.split('+') for x in [produce_names, produce_gate]]
        produce_dict = dict(zip(produce_terms[0][1:], produce_terms[1][1:]))
        return dict(zip([x[0][0] for x in [react_terms, produce_terms]], [react_dict, produce_dict]))

    def get_complexes(self):
        return self.complexes[:]

    def get_top_strands(self):
        return self.top_strands[:]

    def get_top_strand_dict(self):
        t, bm, c = self.params
        in_strands = self.in_strands
        out_dict = self.top_s_dict.copy()
        t_regi = list(range(t+bm-2, t*2 +bm+1))
        out_dict.update(
            dict(
                [ (y, t_regi) for x in self.in_strands for y in x.get_top_strands()]
            ))
        # for in_strand in in_strands:
        #     for seq in in_strand.get_top_strands():
        #         if seq[-2:] in ["-a", "-b"]:
        #             self.top_s_dict.update({seq : list(range(bm+t-2, 1+bm+2*t))})
        #         else:
        #             self.top_s_dict.update({seq : list(range(bm+t-2, 1+bm+2*t))})
        return out_dict

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
        return '{} {} = {}{}: {}\n'.format("component", rxn, comp, param_str, eq)

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
        nr = sum(rxn['stoich_r'])
        if nr not in [0, 1, 2] :
            print("Unexpected stoichiometry!")
            raise Exception("Incorrect reactant stochiometry {} for reaction {}. Must be 0, 1, 2".format(nr, I))
        for reactant_idx in range(nr,2):
            nr += 1
            new_species = 'Fuel' + I
            if reactant_idx == 0:
                new_species = new_species + 'a'
            else:
                new_species = new_species + 'b'
            rxn['stoich_r'].append(1)
            rxn['reactants'].append(new_species)
        reactants = []
        for j, r_species in enumerate(rxn['reactants']):
            if r_species not in species_dict:
                strand = SignalStrand(r_species)
                strand.set_identity_domains(j, rxn_name)
                species_dict[r_species] = strand
            else:
                strand = species_dict[r_species]
            for q in range(rxn['stoich_r'][j]):
                reactants.append(strand)

        np = sum(rxn['stoich_p'])
        if np not in [0, 1, 2] :
            print("Unexpected stoichiometry!")
            raise Exception("Incorrect product stochiometry {} for reaction {}. Must be 0, 1, 2".format(np, I))
        for product_idx in range(np,2):
            np += 1
            new_species = 'Fuel' + I
            if product_idx == 0:
                new_species = new_species + 'c'
            else:
                new_species = new_species + 'd'
            rxn['stoich_p'].append(1)
            rxn['products'].append(new_species)
        products = []
        for j, p_species in enumerate(rxn['products']):
            if p_species not in species_dict:
                strand = SignalStrand(p_species)
                strand.set_identity_domains(j+2, rxn_name)
                species_dict[p_species] = strand
            else:
                strand = species_dict[p_species]
            products.append(strand)
            strand.add_instance(j+2, rxn_name)
            if rxn['stoich_p'][j] == 2:
                products.append(strand)
                strand.add_instance(3, rxn_name)
        gates.append(Bimrxn(rxn_name, reactants, products, d_params))
    strands = list(species_dict.values())
    return (gates, strands)
