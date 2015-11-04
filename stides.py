import numpy as np
import itertools
import logging

class endarray(np.ndarray):
    """
    This is a class for arrays full ends (of type
    adjacent+end+wc-adjacent-of-complementary-end).

    At present, it also handles adjacent+end style ends, but self.end and
    self.comp will return bogus information. It eventually needs to be split up
    into two classes in order to deal with this problem.
    """
    def __new__( cls, array, endtype ):
        if type(array[0]) is str:
            array = np.array([[ nt[x] for x in y ] for y in array ],dtype=np.uint8)
        obj = np.asarray(array,dtype=np.uint8).view(cls)
        obj.endtype = endtype
        return obj
    def __array_finalize__(self, obj):
        if obj is None: return
        self.endtype = getattr(obj, 'endtype', None)
    def _get_ends( self ):
        if self.endtype == 'DT':
            return self[:,:-1]
        elif self.endtype == 'TD':
            return self[:,1:]
    def _get_comps( self ):
        if self.endtype == 'DT':
            return (3-self)[:,::-1][:,:-1]
        elif self.endtype == 'TD':
            return (3-self)[:,::-1][:,1:]
    def _get_adjs( self ):
        if self.endtype == 'DT':
            return self[:,0]
        elif self.endtype == 'TD':
            return self[:,-1]
    def _get_cadjs( self ):
        if self.endtype == 'DT':
            return (3-self)[:,-1]
        elif self.endtype == 'TD':
            return (3-self)[:,0]
    def __len__( self ):
        return self.shape[0]
    def _get_endlen( self ):
        return self.shape[1]-1
    def _get_seqlen( self ):
        return self.shape[1]
    def append(s1,s2):
        assert s1.endtype == s2.endtype
        n = np.vstack( (s1, s2) ).view(endarray)
        n.endtype = s1.endtype
        return n
    def __repr__( self ):
        return "<endarray ({2}): type {0}; {1}>".format( \
                self.endtype, repr(self.tolist()), len(self) )
    def tolist( self ):
        st = ["a","c","g","t"]
        return [ "".join([ st[x] for x in y]) for y in self ]
    endlen = property(_get_endlen)
    seqlen = property(_get_seqlen) # Added for new code compatibility
    ends = property(_get_ends)
    comps = property(_get_comps)
    adjs = property(_get_adjs)
    cadjs = property(_get_cadjs)
    strings = property(tolist)

nt = { 'a': 0, 'c': 1, 'g': 2, 't': 3 }
lton = { 'a': [0],
         'b': [1, 2, 3],
         'c': [1],
         'd': [0, 2, 3],
         'g': [2],
         'h': [0, 1, 3],
         'k': [2, 3],
         'm': [0, 1],
         'n': [0, 1, 2, 3],
         's': [1, 2],
         't': [3],
         'v': [0, 1, 2],
         'w': [0, 3] }
wc = {   'a': 't',
         'b': 'v',
         'c': 'g',
         'd': 'h',
         'g': 'c',
         'h': 'd',
         'k': 'm',
         'm': 'k',
         'n': 'n',
         's': 's',
         't': 'a',
         'v': 'b',
         'w': 'w' }

def values_chunked(items, endtype, chunk_dim=10):
    """
    Given a list of lists of acceptable numbers for each position in a row of
    an array, create every possible row, and return an iterator that returns
    chunks of every possible row up to chunk_dim, iterating dimensions higher
    than chunk_dim.

    Return this as an endarray, with set endtype. This can be easily removed
    for use elsewhere.
    """
    ilengths = [len(x) for x in items]
    n = len(items)
    items = [ np.array(x) for x in items ]
    if n > chunk_dim:
        p = n - chunk_dim
        q = chunk_dim
        outer = itertools.product(*(items[0:p]))
    else:
        p = 0
        q = n
        def outer_iter():
            yield()
        outer = outer_iter()

    chunk = np.zeros([np.prod( ilengths[p:] ),len(items)],dtype=int).view(endarray)
    chunk.endtype = endtype
    chunk[:,p:] = np.indices(ilengths[p:]).reshape(q,-1).T
    for i in range(p,n):
        chunk[:,i] = items[i][chunk[:,i]]
    for seq in outer:
        chunk[:,:p] = seq
        yield chunk

def get_accept_set( endtype, length, interaction, fdev, maxendspurious, spacefilter=None, adjacents=['n','n'], alphabet='n', energyfuncs=None ):
    if not energyfuncs:
        energyfuncs = energyfuncs_santalucia(mismatchtype='max')
    if not spacefilter:
        spacefilter = spacefilter_standard(interaction, interaction*fdev, maxendspurious)
    # Generate the template.
    if endtype == 'DT':
        template = [lton[adjacents[0]]] + [lton[alphabet.lower()]] * length \
                   + [lton[wc[adjacents[1]]]]
    elif endtype == 'TD':
        template = [lton[wc[adjacents[1]]]] + [lton[alphabet.lower()]] * length \
                   + [lton[adjacents[0]]]

    logging.info( "Length {0}, type {1}, adjacents {2}, alphabet {3}.".format( \
                  length,endtype,adjacents,alphabet) )
    logging.debug( "Have template {0}.".format(template,endtype))

    # Create the chunk iterator
    endchunk = values_chunked(template, endtype)

    # Use spacefilter to filter chunks down to usable sequences
    matcharrays = []
    chunknum = 0
    totchunks = None
    totends = np.product( [len(x) for x in template] )
    logging.info( "Have {0} ends in total before any filtering.".format(totends) )
    for chunk in endchunk:
        matcharrays.append(spacefilter(chunk,energyfuncs))
        if not totchunks:
            totchunks = totends/len(chunk)
        chunknum += 1
        logging.debug( "Found {0} filtered ends in chunk {1} of {2}.".format(
            len(matcharrays[-1]), chunknum, totchunks))
    logging.info( "Done with spacefiltering." )
    availends = np.vstack( matcharrays ).view(endarray)
    availends.endtype = endtype
    return availends

def find_end_set_uniform( endtype, length, spacefilter, endfilter, endchooser,\
                          energyfuncs, adjacents=['n','n'], num=0, numtries=1,\
                          oldendfilter=None, oldends = [], alphabet='n'):
    """
        Find a set of ends of uniform length and type satisfying uniform
        constraint functions (eg, constrant functions are the same for each
        end).

        * endtype: right now 'DT' for 3'-terminal ends, and 'TD' for
          5'-terminal ends,
        * length: length of ends, not including adjacents.
        * adjacents (defaults to ['n','n']): acceptable bases for adjacents
          (eg, ['n','n'] or ['c', 'c']) for the ends and their complements,
        * num (defaults to 0): number of ends to find (0 keeps finding until
          available ends are exhausted)
        * numtries (defaults to 1): if > 1, the function will return a list of
          sets of ends that all individually satisfy the constraints, so that
          the best one can be selected manually
        * spacefilter: a "spacefilter" function that takes endarrays and
          filters them down to ends that, not considering spurious
          interactions, are acceptable.
        * endfilter: an "endfilter" function that takes current ends in the
          set, available ends (filtered with current ends), and new ends added,
          and filters the available ends, considering interactions between ends
          (eg, spurious interactions).
        * endchooser: an "endchooser" function that takes current ends in the
          set and available ends, and returns a new end to add to the set.
        * energyfuncs: an "energyfunctions" class that provides the energy
          functions for everything to use.
        * oldends: an endarray of old ends to consider as part of the set
        * alphabet: a single letter specifying what the alphabet for the ends
          should be (eg, four or three-letter code)
        * oldendfilter: a different "endfilter" function for use when filtering
          the available ends using interactions with old ends. This is normally
          not useful, but can be useful if you want, for example, to create a
          sets with higher cross-interactions between two subsets than within
          the two subsets.

        This function is intended to be complicated and featureful. If you want
        something simpler, try easy_ends
    """
    # Generate the template.
    if endtype == 'DT':
        template = [lton[adjacents[0]]] + [lton[alphabet.lower()]] * length \
                   + [lton[wc[adjacents[1]]]]
    elif endtype == 'TD':
        template = [lton[wc[adjacents[1]]]] + [lton[alphabet.lower()]] * length \
                   + [lton[adjacents[0]]]

    logging.info( "Length {0}, type {1}, adjacents {2}, alphabet {3}.".format( \
                  length,endtype,adjacents,alphabet) )
    logging.debug( "Have template {0}.".format(template,endtype))

    # Create the chunk iterator
    endchunk = values_chunked(template, endtype)

    # Use spacefilter to filter chunks down to usable sequences
    matcharrays = []
    chunknum = 0
    totchunks = None
    totends = np.product( [len(x) for x in template] )
    logging.info( "Have {0} ends in total before any filtering.".format(totends) )
    for chunk in endchunk:
        matcharrays.append(spacefilter(chunk,energyfuncs))
        if not totchunks:
            totchunks = totends/len(chunk)
        chunknum += 1
        logging.debug( "Found {0} filtered ends in chunk {1} of {2}.".format(
            len(matcharrays[-1]), chunknum, totchunks))
    logging.info( "Done with spacefiltering." )
    availends = np.vstack( matcharrays ).view(endarray)
    availends.endtype = endtype


    # Use endfilter to filter available sequences taking into account old sequences.
    if len(oldends)>0:
        if type(oldends[0]) is str:
            oldends = endarray(oldends,endtype)
        if oldendfilter:
            availends = oldendfilter( oldends, None, availends, energyfuncs )
        else:
            availends = endfilter( oldends, None, availends, energyfuncs )

    startavail = availends
    endsets = []
    logging.info("Starting with {0} ends.".format(len(availends)))
    while len(endsets) < numtries:
        curends = oldends
        availends = startavail
        numends = 0
        while True:
            newend = endarray( \
                    np.array([endchooser( curends, availends, energyfuncs )]),\
                    endtype)
            logging.debug("Chose end {0}.".format(repr(newend)))
            newend.endtype = endtype
            availends = endfilter( newend, curends, availends, energyfuncs )
            logging.debug("Done filtering.")
            if curends == None:
                curends = newend
            elif len(curends) == 0:
                curends = newend
            else:
                curends = curends.append(newend)
            numends += 1
            logging.info("Now have {0} ends in set, and {1} ends available.".format(\
                    numends,len(availends)))
            if len(availends) == 0:
                logging.warning("Found {0} ends.".format(numends))
                break
            if numends >= num and num > 0:
                break
        endsets.append(curends)

    if len(endsets)>1:
        return endsets
    else:
        return endsets[0]
    # TODO: multiple sets generated, decide which is best.

def enhist( endtype, length, adjacents=['n','n'], alphabet='n',\
            bins=None, energyfuncs=None, plot=False, color='b' ):
    if endtype == 'DT':
        template = [lton[adjacents[0]]] +\
                   [lton[alphabet.lower()]] * length + [lton[wc[adjacents[1]]]]
    elif endtype == 'TD':
        template = [lton[wc[adjacents[1]]]] +\
                   [lton[alphabet.lower()]] * length + [lton[adjacents[0]]]

    if not energyfuncs:
        energyfuncs = energyfuncs_santalucia(mismatchtype='max')

    minbin = 0.8*energyfuncs.matching_uniform( \
            endarray( [([0,3]*length)[0:length+2]], endtype ) )
    maxbin = 1.1*energyfuncs.matching_uniform( \
            endarray( [([1,2]*length)[0:length+2]], endtype ) )

    if not bins:
        bins = np.arange(minbin,maxbin,0.1)

    logging.debug( "Have template {0} and type {1}.".format(template,endtype))

    # Create the chunk iterator
    endchunk = values_chunked(template, endtype)

    hist = np.zeros(len(bins)-1,dtype='int')
    totends = np.product( [len(x) for x in template] )
    finishedends = 0
    info = {'min': np.inf, 'max': -np.inf, 'mean': 0}
    for chunk in endchunk:
        matchens = energyfuncs.matching_uniform(chunk)
        hist += np.histogram(matchens,bins)[0]
        info['max'] = max( info['max'], np.amax(matchens) )
        info['min'] = min( info['min'], np.amin(matchens) )
        info['mean'] = info['mean']*(finishedends)/(len(chunk)+finishedends)\
                + np.mean(matchens)*len(chunk)/(len(chunk)+finishedends)
        finishedends += len(matchens)
        logging.debug( "Done with {0}/{1} ends.".format(finishedends,totends) )

    x = (bins[:-1]+bins[1:])/2
    n = hist
    info['emean'] = np.sum( n*x, dtype='double' ) / np.sum( n, dtype='int64' )
    info['estd'] =  np.sqrt( np.sum( n*(x-info['emean'])**2, dtype='double' )\
            / np.sum( n, dtype='int64' ) )
    cs = np.cumsum(n)
    info['emedian'] = x[ np.flatnonzero( cs >= cs[-1]/2.0 )[0] ]

    if plot:
        import matplotlib.pyplot as plt

        plt.bar( bins[:-1], hist, width=(bins[1]-bins[0]),\
                label="Type {3}, Length {0}, Adjs {1}, Alph {2}".format(\
                length,adjacents,alphabet,endtype),color=color )
        plt.title(\
            "Matching Energies of Ends of Type {3}, Length {0}, Adjs {1}, Alph {2}".format(\
            length,adjacents,alphabet,endtype))
        plt.xlabel("Standard Free Energy (-kcal/mol)")
        plt.ylabel("Number of Ends")
        #plt.show()

    return (hist,bins,info)


def easyends( endtype, endlength, number=0, interaction=None, fdev=0.05,\
              maxspurious=0.5, maxendspurious=None, tries=1, oldends=[],\
              adjs=['n','n'], energyfuncs=None, alphabet='n', echoose=None ):
    """
    Easyends is an attempt at creating an easy-to-use function for finding sets
    of ends.

    * endtype: specifies the type of end being considered. The system for
      classifying end types goes from 5' to 3', and consists of letters
      describing each side of the end. For example, an end that starts after a
      double-stranded region on the 5' side and ends at the end of the strand
      would be 'DT', while one that starts at the beginning of a strand on the
      5' side and ends in a double-stranded region would be 'TD'. 'T' stands
      for terminal, 'D' stands for double-stranded region, and 'S' stands for
      single-stranded region. 'S', however, is not currently supported.
    * endlength: specifies the length of end being considered, not including
      adjacent bases.
    * number (optional): specifies the number of ends to find.  If zero or not
      provided, easyends tries to find as many ends as possible.
    * interaction (optional): a positive number corresponding to the desired
      standard free energy for hybridization of matching sticky ends. If not
      provided, easyends calculates an optimal value based on the sequence
      space.
    * fdev (default 0.05): the fractional deviation (above or below) of
      allowable matching energies.  maxspurious (default 0.5): the maximum
      spurious interaction, as a fraction of the matching interaction.
    * maxendspurious (default None): if provided, maxspurious is only used for
      spurious interactions between ends defined as ends, and ends defined as
      complements. Maxendspurious is then the maximum spurious interaction
      between ends and ends, and complements and complements. In a system
      where spurious interactions between ends and complements are more important
      than other spurious interactions, this can allow for better sets of ends.
    * tries (default 1): if > 1, easyends will return a list of sets of ends,
      all satisfying the constraints.
    * oldends (optional): a list of ends to be considered as already part of
      the set.
    * adjacents (default ['n','n']): allowable adjacent bases for ends and
      complements.
    * energyfuncs (optional): an energyfuncs class providing the energy
      calculation functions. You probably don't need to change this.
    * alphabet (default 'n'): The alphabet to use for ends, allowing
      for three-letter codes.
    """

    if not energyfuncs:
        efunc = energyfuncs_santalucia(mismatchtype='max')
    else:
        efunc=energyfuncs
    if (not interaction) or (interaction == 0):
        interaction = enhist( endtype, endlength, energyfuncs=efunc,\
                adjacents=adjs, alphabet=alphabet )[2]['emedian']
        logging.warning("Calculated optimal interaction energy is {0}.".format(interaction))
    maxcompspurious = maxspurious*interaction
    if not maxendspurious:
        maxendspurious = maxspurious*interaction
    else:
        maxendspurious = maxendspurious*interaction

    sfilt = spacefilter_standard(interaction, interaction*fdev, maxendspurious)
    efilt = endfilter_standard_advanced(maxcompspurious,maxendspurious)
    if not echoose:
        echoose = endchooser_standard(interaction)

    return find_end_set_uniform( endtype, endlength, sfilt, efilt, echoose,\
           energyfuncs=efunc, numtries=tries, oldends=oldends,adjacents=adjs,\
           num=number,alphabet=alphabet)

def easy_space( endtype, endlength, interaction=None, fdev=0.05,\
              maxspurious=0.5, maxendspurious=None, tries=1, oldends=[],\
              adjs=['n','n'], energyfuncs=None, alphabet='n', echoose=None, runnx=False ):
    length = endlength
    if not energyfuncs:
        efunc = energyfuncs_santalucia(mismatchtype='max')
        energyfuncs = efunc
    else:
        efunc=energyfuncs
    if (not interaction) or (interaction == 0):
        interaction = enhist( endtype, endlength, energyfuncs=efunc,\
                adjacents=adjs, alphabet=alphabet )[2]['emedian']
        logging.warning("Calculated optimal interaction energy is {0}.".format(interaction))
    maxcompspurious = maxspurious*interaction
    if not maxendspurious:
        maxendspurious = maxspurious*interaction
    else:
        maxendspurious = maxendspurious*interaction

    sfilt = spacefilter_standard(interaction, interaction*fdev, maxendspurious)
    spacefilter = sfilt
    efilt = endfilter_standard_advanced(maxcompspurious,maxendspurious)
    if not echoose:
        echoose = endchooser_standard(interaction)

    adjacents = adjs

    if endtype == 'DT':
        template = [lton[adjacents[0]]] + [lton[alphabet.lower()]] * length \
                   + [lton[wc[adjacents[1]]]]
    elif endtype == 'TD':
        template = [lton[wc[adjacents[1]]]] + [lton[alphabet.lower()]] * length \
                   + [lton[adjacents[0]]]

    # Create the chunk iterator
    endchunk = values_chunked(template, endtype)

    # Use spacefilter to filter chunks down to usable sequences
    matcharrays = []
    chunknum = 0
    totchunks = None
    totends = np.product( [len(x) for x in template] )
    logging.info( "Have {0} ends in total before any filtering.".format(totends) )
    for chunk in endchunk:
        matcharrays.append(spacefilter(chunk,energyfuncs))
        if not totchunks:
            totchunks = totends/len(chunk)
        chunknum += 1
        logging.debug( "Found {0} filtered ends in chunk {1} of {2}.".format(
            len(matcharrays[-1]), chunknum, totchunks))
    logging.info( "Done with spacefiltering." )
    availends = np.vstack( matcharrays ).view(endarray)
    availends.endtype = endtype

    availendsr = np.repeat(availends,len(availends), axis=0)
    availendst = np.tile(availends,(len(availends),1))

    vals_ee = energyfuncs.uniform( availendsr.ends, availendst.ends )
    vals_ec = energyfuncs.uniform( availendsr.ends, availendst.comps )
    vals_ce = energyfuncs.uniform( availendsr.comps, availendst.ends )
    vals_cc = energyfuncs.uniform( availendsr.comps, availendst.comps )

    vals_tf = ( (vals_ee < maxendspurious) & (vals_cc < maxendspurious) & (vals_ec < maxcompspurious) & (vals_ce < maxcompspurious) )
    zipendsnf = zip(availendsr.tolist(),availendst.tolist())
    zipends = [ zipendsnf[x] for x in np.flatnonzero(vals_tf) ]

    if not runnx:
        return zipends

    import networkx as nx
    G = nx.Graph()
    G.add_edges_from(zipends)

    maxl = 0
    for clique in nx.clique.find_cliques(G):
        if len(clique) > maxl:
            print "Found clique length {0}: {1}".format(len(clique),clique)
            maxl = len(clique)
            maxc = clique
    return maxc

def spacefilter_standard(desint, dev, maxself):
    """
    A spacefilter function: filters to ends that have a end-complement
    interaction of between desint-dev and desint+dev, and a self-interaction
    (end-end or comp-comp) of less than maxself.
    """
    def spacefilter( fullends, energyfuncs ):
        matchenergies = energyfuncs.matching_uniform( fullends )
        i = np.flatnonzero( (matchenergies<desint+dev) & (matchenergies>desint-dev) )
        matchenergies = matchenergies[i]; fullends = fullends[i]
        selfselfenergies = energyfuncs.uniform( fullends.ends, fullends.ends )
        compcompenergies = energyfuncs.uniform( fullends.comps, fullends.comps )
        i = np.flatnonzero( (selfselfenergies < maxself) & (compcompenergies < maxself) )
        return fullends[i]
    return spacefilter

def endfilter_standard(maxspurious):
    """
    An endfilter function: filters out ends that have any (end-end, end-comp,
    comp-end, comp-comp) interactions with new ends above maxspurious.
    """
    def endfilter( newends, currentends, availends, energyfuncs ):
        endendspurious = energyfuncs.uniform(
                np.repeat(newends.ends,len(availends),0),
                np.tile(availends.ends,(len(newends),1)) ).reshape(
                        len(availends), len(newends), order='F' )
        endcompspurious = energyfuncs.uniform(
                np.repeat(newends.ends,len(availends),0),
                np.tile(availends.comps,(len(newends),1)) ).reshape(
                        len(availends), len(newends), order='F')
        compendspurious = energyfuncs.uniform(
                np.repeat(newends.comps,len(availends),0),
                np.tile(availends.ends,(len(newends),1)) ).reshape(
                        len(availends), len(newends), order='F' )
        compcompspurious = energyfuncs.uniform(
                np.repeat(newends.comps,len(availends),0),
                np.tile(availends.comps,(len(newends),1)) ).reshape(
                        len(availends), len(newends), order='F' )
        highspurious = np.amax( np.hstack( (endendspurious, compendspurious,
            endcompspurious, compcompspurious) ), 1 )
        #pdb.set_trace()
        return availends[ highspurious < maxspurious ]
    return endfilter

def endfilter_standard_advanced(maxcompspurious,maxendspurious):
    """
    An endfilter function: filters out ends that have end-comp or comp-end
    interactions above maxcompspurious, and end-end or comp-comp interactions
    above maxendspurious.
    """
    def endfilter( newends, currentends, availends, energyfuncs ):
        endendspurious = energyfuncs.uniform(
                np.repeat(newends.ends,len(availends),0),
                np.tile(availends.ends,(len(newends),1)) ).reshape(
                        len(availends), len(newends), order='F' )
        endcompspurious = energyfuncs.uniform(
                np.repeat(newends.ends,len(availends),0),
                np.tile(availends.comps,(len(newends),1)) ).reshape(
                        len(availends), len(newends), order='F')
        compendspurious = energyfuncs.uniform(
                np.repeat(newends.comps,len(availends),0),
                np.tile(availends.ends,(len(newends),1)) ).reshape(
                        len(availends), len(newends), order='F' )
        compcompspurious = energyfuncs.uniform(
                np.repeat(newends.comps,len(availends),0),
                np.tile(availends.comps,(len(newends),1)) ).reshape(
                        len(availends), len(newends), order='F' )
        highendspurious = np.amax( np.hstack( (endendspurious,
            compcompspurious) ), 1 )
        highcompspurious = np.amax( np.hstack( (compendspurious,
            endcompspurious) ), 1 )
        #pdb.set_trace()
        return availends[ (highendspurious < maxendspurious) &
                (highcompspurious < maxcompspurious) ]
    return endfilter

def energy_array_uniform( seqs, energyfuncs ):
    """
    Given an endarray and a set of sequences, return an array of the
    interactions between them, including their complements.
    """
    seqs = seqs.ends.append(seqs.comps)
    return energyfuncs.uniform( np.repeat(seqs,seqs.shape[0],0),
            np.tile(seqs,(seqs.shape[0],1)) ).reshape(
                    (seqs.shape[0],seqs.shape[0]) )

def endchooser_standard(desint):
    """
    An endchooser function: return a random end with end-comp energy closest to
    desint.
    """
    def endchooser( currentends, availends, energyfuncs ):
        ddiff = ( energyfuncs.matching_uniform( availends ) - desint )**2
        choices = np.flatnonzero( ddiff == np.amin(ddiff) )
        newend = availends[choices[np.random.randint(0,len(choices))]]
        return newend
    return endchooser

def endchooser_random():
    """
    An endchooser function: return a random end with end-comp energy closest to
    desint.
    """
    def endchooser( currentends, availends, energyfuncs ):
        newend = availends[np.random.randint(0,len(availends))]
        return newend
    return endchooser

class energyfuncs_santalucia:
    """
    Energy functions based on SantaLucia's 2004 paper.

    mismatchtype is one of 'max', 'loop', or 'dangle', specifying how to
    consider mismatches.  'max' is probably the best choice, but is slowest -
    it takes the maximum interaction of the 'loop' and 'dangle' options.
    """
    def __init__(self, mismatchtype='max'):
        import os
        try:
            import pkg_resources
            dsb = pkg_resources.resource_stream(__name__, os.path.join('params','dnastackingbig.csv'))
        except:
            try:
                this_dir, this_filename = os.path.split(__file__)
                dsb = open( os.path.join(this_dir, "params", "dnastackingbig.csv") )
            except IOError:
                raise IOError("Error loading dnastackingbig.csv")
        self.nndG_full = -np.loadtxt(dsb ,delimiter=',')
        dsb.close()
        self.initdG = 0.0 # 1.96 DISABLED FOR NOW
        self.nndG = self.nndG_full[np.arange(0,16),15-np.arange(0,16)]
        if mismatchtype == 'max':
            self.uniform = lambda x,y: np.maximum( self.uniform_loopmismatch(x,y), \
                                                   self.uniform_danglemismatch(x,y) \
                                                 )
        elif mismatchtype == 'loop':
            self.uniform = self.uniform_loopmismatch
        elif mismatchtype == 'dangle':
            self.uniform = self.uniform_danglemismatch
        else:
            raise InputError("Mismatchtype {0} is not supported.".format(mismatchtype))

    def matching_uniform(self, seqs):
        return np.sum(self.nndG[tops(seqs)],axis=1) - self.initdG

    def uniform_loopmismatch(self, seqs1, seqs2):
        if seqs1.shape != seqs2.shape:
            if seqs1.ndim == 1:
                seqs1 = endarray( np.repeat(np.array([seqs1]),seqs2.shape[0],0), seqs1.endtype )
            else:
                raise InputError("Lengths of sequence arrays are not acceptable.")
        assert seqs1.endtype == seqs2.endtype
        endtype = seqs1.endtype

        endlen = seqs1.endlen
        plen = endlen-1

        # TODO: replace this with cleaner code
        if endtype=='DT':
            ps1 = seqs1[:,1:-1]*4+seqs1[:,2:]
            pa1 = seqs1[:,0]*4+seqs1[:,1]
            pac1 = (3-seqs1[:,0])*4+seqs2[:,-1]
            ps2 = seqs2[:,::-1][:,:-2]*4+seqs2[:,::-1][:,1:-1]
            pa2 = seqs2[:,0]*4+seqs2[:,1]
            pac2 = (3-seqs2[:,0])*4+seqs1[:,-1]
        if endtype=='TD':
            ps1 = seqs1[:,:-2]*4+seqs1[:,1:-1]
            pa1 = seqs1[:,-2]*4+seqs1[:,-1]
            pac1 = seqs2[:,0]*4+(3-seqs1[:,-1])
            ps2 = seqs2[:,::-1][:,1:-1]*4+seqs2[:,::-1][:,2:]
            pa2 = seqs2[:,-2]*4+seqs2[:,-1]
            pac2 = (seqs1[:,-1])*4+(3-seqs2[:,-1])

        # Shift here is considering the first strand as fixed, and the second one as
        # shifting.  The shift is the offset of the bottom one in terms of pair
        # sequences (thus +2 and -1 instead of +1 and 0).
        en = np.zeros( (ps1.shape[0], 2*plen) )
        for shift in range(-plen+1,plen):
            #import pdb
            #pdb.set_trace()
            en[:,plen+shift-1] = np.sum( \
                    self.nndG_full[ ps1[:,max(shift,0):plen+shift], \
                               ps2[:,max(-shift,0):plen-shift] ], \
                               axis=1)
        en[:,plen-1] = en[:,plen-1] + self.nndG_full[pa1,pac1] + self.nndG_full[pa2,pac2]
        return np.amax(en,1) - self.initdG

    def uniform_danglemismatch(self, seqs1,seqs2,fast=True):
        if seqs1.shape != seqs2.shape:
            if seqs1.ndim == 1:
                seqs1 = endarray( np.repeat(np.array([seqs1]),seqs2.shape[0],0), seqs1.endtype )
            else:
                raise InputError("Lengths of sequence arrays are not acceptable.")
        assert seqs1.endtype == seqs2.endtype
        endtype = seqs1.endtype
        s1 = tops(seqs1)
        s2 = tops(seqs2)
        l = s1.shape[1]
        s2r = np.fliplr(np.invert(s2)%16)
        s2r = s2r/4 + 4*(s2r%4)
        m = np.zeros((s1.shape[0],2*np.sum(np.arange(2,l+1))+l+1))
        r = np.zeros(m.shape[0])
        z = 0;
        if endtype == 'TD':
            s1c = s1[:,0:-1]
            s2rc = s2r[:,1:]
            s1l = np.hstack(( (4*(s2r[:,0]/4) + s1[:,0]/4).reshape(-1,1) , s1 ))
            s2rl = np.hstack(( s2r , (4*(s2r[:,-1]%4) + s1[:,-1]%4).reshape(-1,1) ))
        elif endtype == 'DT':
            s1c = s1[:,1:]
            s2rc = s2r[:,0:-1]
            s2rl = np.hstack(( (4*(s1[:,0]/4) + s2r[:,0]/4).reshape(-1,1) , s2r ))
            s1l = np.hstack(( s1 , (4*(s1[:,-1]%4) + s2r[:,-1]%4).reshape(-1,1) ))
        for o in range(1,l-1):
            zn = l-1-o
            m[:,z:z+zn] = ( s1c[:,:-o]==s2rc[:,o:] ) * self.nndG[s1c[:,:-o]]
            z = z+zn+2
            m[:,z:z+zn] = ( s2rc[:,:-o]==s1c[:,o:] ) * self.nndG[s2rc[:,:-o]]
            z = z+zn+2
        m[:,z:z+l+1] = (s1l == s2rl) * self.nndG[s1l]
        i = 0
        im = len(m)
        # This needs to be changed to something faster
        if not fast:
            for xi in range(0,m.shape[0]):
                gm = 0
                g = 0
                for y in m[xi,:]:
                    if y == 0:
                        g = 0
                    else:
                        g += y
                        if gm > g:
                            gm = g
                r[xi] = gm
                i+=1
                if not i%1000:
                    print "%d/%d" % (i,im)
        else:
            import _stickyext
            x = m
            _stickyext.fastsub(x,r)

        return r-self.initdG


tops = lambda s: 4*s[:,:-1]+s[:,1:]


if __name__ == '__main__':
    # OK, we're being run as a script, so let's deal with it.

    # Right now this is very rudimentary

    import sys
    args = [ eval(x) for x in sys.argv[1:] ]

    E = easyends( *args )

    print repr(E)
