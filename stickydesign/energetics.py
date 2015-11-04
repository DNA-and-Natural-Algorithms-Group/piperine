from .stickydesign import *

def tops(s):
    return 4*s[:,:-1]+s[:,1:]

class energetics_santalucia:
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


