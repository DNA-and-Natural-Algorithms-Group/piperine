from __future__ import division, print_function
import sys
from pkg_resources import Requirement, resource_stream
import numpy as np
import itertools
import logging

from . import translation

nt = { 'a': 0, 'c': 1, 'g': 2, 't': 3 }
tops = lambda s: 4*s[:,:-1]+s[:,1:]

class ToeholdSpecificationError(ValueError):
    '''
    Error raised when stickydesign cannot satisfy user's toehold specifications.
    '''
    def __init__(self, message):
        #self.expression = expression
        self.message = message

def exceptionhook(exception_type, exception, traceback, default_hook = sys.excepthook):
    if 'Toehold' in exception_type.__name__ :
        print("{}: {}".format('RuntimeError', exception.message))
    else:
        default_hook(exception_type, exception, traceback)

sys.excepthook = exceptionhook

class energyfuncs:
    """
    Energy functions based on SantaLucia's 2004 paper.

    mismatchtype is one of 'max', 'loop', or 'dangle', specifying how to
    consider mismatches.  'max' is probably the best choice, but is slowest -
    it takes the maximum interaction of the 'loop' and 'dangle' options.
    """
    def __init__(self, targetdG=7, length=7, deviation=0.5, max_spurious=0.4):
        import os
        try:
            dsb = resource_stream('stickydesign', 'stickydesign/params/dnastackingbig.csv')
        except:
            try:
                dsb = resource_stream('stickydesign', 'params/dnastackingbig.csv')
            except IOError:
                raise IOError("Error loading dnastackingbig.csv")
        try:
            dgl = resource_stream('piperine', 'data/dnadangle.csv')
        except:
            try:
                this_dir, this_filename = os.path.split(__file__)
                dgl = open( os.path.join(this_dir, "data", "dnadangle.csv") )
            except IOError:
                raise IOError("Error loading dnadangle.csv")
        self.targetdG=targetdG
        self.deviation=deviation
        self.max_spurious=max_spurious
        self.length=length
        self.nndG_full = -np.loadtxt(dsb ,delimiter=',')
        self.dgldG_full = -np.loadtxt(dgl ,delimiter=',')
        self.taildG = 1.3
        dsb.close()
        dgl.close()
        self.initdG = 0.0 # 1.96 DISABLED FOR NOW
        self.nndG = self.nndG_full[np.arange(0,16),15-np.arange(0,16)]
        # 30-01-15: The only dangle contexts we are interested in are 3' dangle
        # s. Select those from the Santa Lucia table. We'll have to flip the
        # order of the vector, though, to mach the 5->3 orientation of the gene
        # rated toeholds. To flip, we need to count up to 15 with the opposite-
        # endian order, in terms of quaternary representation
        indcs = 4*np.tile(np.arange(4), 4) + np.repeat(np.arange(4), 4)
        self.dgldG = self.dgldG_full[1, indcs]
        # 30-01-15: As of now, the dangle base is set to C, so make a lookup ta
        # ble ordered by terminating toehold base
        self.dgldG_fixedC = self.dgldG_full[1, np.arange(4) + 4 * nt['c']]
        self.uniform = lambda x,y: np.maximum( self.uniform_loopmismatch(x,y), \
                                               self.uniform_danglemismatch(x,y) \
                                             )

    def th_external_3_dG(self, seqs):
        # Convert nearest-neighbor stacks to dG-table lookup indices.
        # Sum up the near-neighbor energy contributions
        # Add context-specific dG values, eg tail or dangle contributions
        seqs_len = np.size(seqs, 1)
        # The external context involves a 3' dangle, so exclude the 3' flank
        # base.
        cols_external = np.arange(seqs_len-1)
        tops_external = tops(seqs[:, cols_external])
        nndG_external = np.sum(self.nndG[tops_external], 1)
        return nndG_external - self.taildG - self.initdG

    def th_external_5_dG(self, seqs):
        # Convert nearest-neighbor stacks to dG-table lookup indices.
        # Sum up the near-neighbor energy contributions
        # Add context-specific dG values, eg tail or dangle contributions
        seqs_len = np.size(seqs, 1)
        # The external context involves a 5' dangle, so exclude the 5' flank
        # base.
        cols_external = np.arange(1, seqs_len)
        tops_external = tops(seqs[:, cols_external])
        nndG_external = np.sum(self.nndG[tops_external], 1)
        return nndG_external - self.taildG - self.initdG

    def th_external_dG(self, seqs):
        # Make a boolean vector representing which toeholds' external 3' context dG
        # is further from the target dG than than their external 5' context dG
        dG_external3 = self.th_external_3_dG(seqs)
        dG_external5 = self.th_external_5_dG(seqs)
        external_further_bool = np.abs(dG_external3 - self.targetdG) >\
                                np.abs(dG_external5 - self.targetdG)
        return np.choose(external_further_bool, [dG_external5, dG_external3])

    def th_internal_dG(self, seqs):
        # Convert nearest-neighbor stacks to dG-table lookup indices
        # Sum up and return the near-neighbor energy contributions
        # Add context-specific dG values, eg tail or dangle contributions
        seqs_len = np.size(seqs, 1)
        tops_internal = tops(seqs)
        nndG_internal = np.sum(self.nndG[tops_internal], 1)
        return nndG_internal - self.taildG - self.initdG

    def matching_uniform(self, seqs):
        # Make a boolean vector representing which toeholds' external context dG
        # is further from the target dG than than their internal context dG
        dG_external = self.th_external_dG(seqs)
        dG_internal = self.th_internal_dG(seqs)
        external_further_bool = np.abs(dG_external - self.targetdG) >\
                                np.abs(dG_internal - self.targetdG)
        return np.choose(external_further_bool, [dG_internal, dG_external])

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

        # Run through the
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
        s2r = s2r//4 + 4*(s2r%4)
        m = np.zeros((s1.shape[0],2*np.sum(np.arange(2,l+1))+l+1))
        r = np.zeros(m.shape[0])
        z = 0;
        if endtype == 'TD':
            s1c = s1[:,0:-1]
            s2rc = s2r[:,1:]
            s1l = np.hstack(( (4*(s2r[:,0]//4) + s1[:,0]//4).reshape(-1,1) , s1 ))
            s2rl = np.hstack(( s2r , (4*(s2r[:,-1]%4) + s1[:,-1]%4).reshape(-1,1) ))
        elif endtype == 'DT':
            s1c = s1[:,1:]
            s2rc = s2r[:,0:-1]
            s2rl = np.hstack(( (4*(s1[:,0]//4) + s2r[:,0]//4).reshape(-1,1) , s2r ))
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
                    print("%d/%d" % (i,im))
        else:
            from stickydesign import _stickyext
            x = m
            _stickyext.fastsub(x,r)

        return r-self.initdG

    def score_toeholds(self, toeholds):
        import stickydesign as sd
        toeholds = translation.flatten(toeholds)
        toeholds_flanked = [ 'c' + th.lower() + 'c' for th in toeholds]
        ends = sd.endarray(toeholds_flanked, 'TD')
        e_fn_list = [self.th_external_3_dG, self.th_external_5_dG,
                     self.th_internal_dG, self.th_external_dG]
        e_vec_list = [ fn(ends) for fn in e_fn_list ]
        e_vec_all = np.concatenate( e_vec_list )
        e_err = np.abs(e_vec_all.mean() - self.targetdG)
        e_rng = e_vec_all.max() - e_vec_all.min()
        return (e_err, e_rng)

    def calculate_unrestricted_toehold_characteristics(self):
        import stickydesign as sd
        ends = sd.easyends('TD',
                           self.length,
                           alphabet='h',
                           adjs=['c', 'g'],
                           energetics=self)
        n_ends = len(ends)
        e_array = sd.energy_array_uniform(ends, self)
        e_array = e_array[n_ends:, :n_ends]
        for i in range(n_ends):
            e_array[i,i] = 0
        e_spr = e_array.max()/self.targetdG
        e_vec_ext = self.th_external_dG(ends)
        e_vec_int = self.th_internal_dG(ends)
        e_vec_all = np.concatenate( (e_vec_int, e_vec_ext))
        e_avg = e_vec_all.mean()
        e_dev = np.max(np.abs(e_vec_all - self.targetdG))
        return e_avg, e_spr, e_dev, n_ends

    def get_toeholds(self, n_ths=6, timeout=2):
        from  time import time
        import stickydesign as sd
        """ Generate specified stickyends for the Soloveichik DSD approach

        A given run of stickydesign may not generate toeholds that match the
        Soloveichik approach. This function reports whether the run was successful
        or not. Otherwise, the toeholds are matched to respect the backwards
        strand back-to-back toeholds.

        Args:
            n_ths: Number of toeholds to generate. (6)
            timeout: Time duration allowed for finding toeholds in seconds. (8)
        Returns:
            List of toehold strings
        """
        # Give StickyDesign a set of trivial, single-nucleotide toeholds to avoid poor
        # designs. I'm not sure if this helps now, but it did once.
        avoid_list = [i * int(self.length + 2) for i in ['a', 't']]

        # Generate toeholds
        fdev = self.deviation / self.targetdG
        notoes = True
        startime = time()
        while notoes:
            try:
                ends = sd.easyends('TD',
                                   self.length,
                                   interaction=self.targetdG,
                                   fdev=fdev,
                                   alphabet='h',
                                   adjs=['c', 'g'],
                                   maxspurious=self.max_spurious,
                                   energetics=self,
                                   oldends=avoid_list)
                notoes = len(ends) < n_ths + len(avoid_list)
                if (time() - startime) > timeout:
                    e_avg, e_spr, e_dev, n_ends = self.calculate_unrestricted_toehold_characteristics()
                    msg = "Cannot make toeholds to user specification! Try target energy:{:.2}, maxspurious:{:.2}, deviation:{:.2}, which makes {:d} toeholds."
                    exception = ToeholdSpecificationError(msg.format(e_avg, e_spr, e_dev, n_ends))
                    raise exception
            except ValueError:
                e_avg, e_spr, e_dev, n_ends = self.calculate_unrestricted_toehold_characteristics()
                msg = "Cannot make toeholds to user specification! Try target energy:{:.2}, maxspurious:{:.2}, deviation:{:.2}, which makes {:d} toeholds."
                exception = ToeholdSpecificationError(msg.format(e_avg, e_spr, e_dev, n_ends))
                raise exception

        th_cands = ends.tolist()
        # remove "avoid" sequences
        th_cands = th_cands[len(avoid_list):]
        # Make as many end in c as possible
        th_cands = th_cands[:n_ths]
        th_all = [ th[1:-1] for th in th_cands]
        ends_full = sd.endarray(th_cands, 'TD')
        ends_all = sd.endarray(th_all, 'TD')
        return ends_all.tolist()


