#!/usr/bin/env python
# Ok, here we go. Whether a library or not, we need to set up ctypes

import ctypes
import numpy as np
import logging as log
import random

group = {"A": "A", "T": "T", "C": "C", "G": "G",
         "W": "AT", "S": "CG", "M": "AC", "K": "GT", 
         "B": "CGT", "V": "ACG", "D": "AGT", "H": "ACT",
         "N": "ACGT", " ": ""}


binbase = { ' ': 0, 'A': 1, 'C': 2, 'G': 3, 'T': 4 }
binbaser = dict( (v,k) for k,v in binbase.items() )
wcletterdict = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
bindict = np.array( [ 0b0000, 0b0001, 0b0010, 0b0100, 0b1000 ]  )
binrevd = dict( (v,k) for k,v in enumerate(bindict) )
letterfrombind = { 0b0000: ' ', 0b0001: 'A', 0b0010: 'C', 0b0100: 'G', 0b1000: 'T'}

def binrev(array):
	return np.array([ binrevd[x] for x in array ])

def letterfrombin(array):
	return np.array([ letterfrombind[x] for x in array ])

def tobin( strarray ):
	binarray = np.zeros(len(strarray),dtype='int')

	for i in range(0,len(strarray)):
		for b in strarray[i]:
			binarray[i] = binarray[i] | bindict[binbase[b]]
	
	return binarray

orange = np.arange(0,16)
wcbin = np.array( [ sum( 1<<(4-1-i) for i in range(4) if o>>i&1 ) for o in orange ] )

ncbin = (1&orange) + (1&orange>>1) + (1&orange>>2) + (1&orange>>3)

choicearray = np.zeros(orange.shape, dtype='object')
choicearraybin = np.zeros(orange.shape, dtype='object')

for x in orange:
	choicearray[x] = set()
	choicearraybin[x] = set()
	if (1&orange[x]): choicearray[x].add('A') 
	if (1&orange[x]>>1): choicearray[x].add('C') 
	if (1&orange[x]>>2): choicearray[x].add('G') 
	if (1&orange[x]>>3): choicearray[x].add('T') 
	if (1&orange[x]): choicearraybin[x].add(1) 
	if (1&orange[x]>>1): choicearraybin[x].add(2) 
	if (1&orange[x]>>2): choicearraybin[x].add(3) 
	if (1&orange[x]>>3): choicearraybin[x].add(4) 

class SpuriousScorer:
	def __init__(self, wc, eq):
		self.N = len(eq)
		self.reltype = self.N * ctypes.c_int
		self.eq_aref = ctypes.byref(self.reltype( *eq.astype('int') ))
		self.wc_aref = ctypes.byref(self.reltype( *wc.astype('int') ))
		self.spuriouslib = ctypes.CDLL('spuriousC.so')
		self._spurious = self.spuriouslib.spurious

	def spurious(self, seq, kmin=3, kmax=8, mismatch=0):
		restype = ctypes.POINTER( (kmax-kmin+1)*3*ctypes.c_int )
		self._spurious.restype = restype
		c = self._spurious( mismatch, seq, kmin, kmax, self.wc_aref, self.eq_aref )
		return np.array([x for x in c[0]]).reshape((3,-1))

def bitand(array):
	z = 15
	for x in array:
		z = z & x
	return z

class SpuriousSystemNum:
	def spurious(self, seq, kmin=3, kmax=8, mismatch=0):
		restype = ctypes.POINTER( (kmax-kmin+1)*3*ctypes.c_int )
		self._spurious.restype = restype
		c = self._spurious( mismatch, seq.tostring(), self.N, kmin, kmax, self.wc_aref, self.eq_aref )
		return np.array([x for x in c[0]]).reshape((3,-1))

	def __init__(self, templ, wc, eq, beta=5.0, spurious_weights=np.array([[1.0,1.0,1.0]])):

		# first things first, let's find out what our free parameters
		# unconstrained parameters have eq[i] == i+1 (equal to itself), and wc[i] > i (wc with something later, not earlier)
		self.N = len(wc)
		self.reltype = self.N * ctypes.c_int
		self.eq_aref = ctypes.byref(self.reltype( *eq.astype('int') ))
		self.wc_aref = ctypes.byref(self.reltype( *wc.astype('int') ))
		self.spuriouslib = ctypes.CDLL('./spuriousC.so')
		self._spurious = self.spuriouslib.spuriousnum
		self.spuriouslib.make_bad6mer()
		self.spuriouslib.score_verboten_num.restype = ctypes.c_double
		self.index = np.arange(1,self.N+1)
		self.unconstrained = (eq == self.index) & ((wc > self.index) | (wc == -1))
		self._beta = beta
		self._spurious_weights = spurious_weights
		self.weightarray = np.outer(spurious_weights,beta**np.array([0,1,2,3,4,5]))

		self.eq = eq.astype('int')
		self.wc = wc.astype('int')

		# now let's convert our template
		self.choices = tobin(np.array([ group[x] for x in templ ]))

		for (i,ii) in zip(np.arange(0,self.N),self.index):

			self.choices[i] = self.choices[i] & bitand(self.choices[self.eq == ii])
			self.choices[i] = self.choices[i] & bitand(wcbin[self.choices[self.wc == ii]].copy())				

		self.nchoices = ncbin[self.choices].copy()

		# Now, for items that don't have choices, we need to fix them and anything that depends on them.

		#seq[(self.nchoices==1)] = binrev(self.choices[(self.nchoices==1)])

		# actually free parameters have choices and are unconstrained
		self.free = (self.unconstrained) & (self.nchoices>1)
		self.blank = self.nchoices==0
		self.constrained = (~self.unconstrained) & (~self.blank)
		self.set = (self.unconstrained) & (self.nchoices==1)
		self.freeindex = self.index[self.free]-1
		self.lenfree = len(self.freeindex)

		log.info("%d unconstrained bases out of %d. %d have choices." % (np.sum(self.unconstrained),self.N,np.sum(self.free)) )

	def constrain( self, seq ):
		"""Fix a new sequence up to fit constraints."""
		seq[self.set] = binrev(self.choices[(self.set)])
		pass

	def spurious_score(self, seq):
		return np.sum( self.spurious(seq) * self.weightarray ) + \
			np.sum( self.spurious(seq,mismatch=1,kmin=5,kmax=10) * self.weightarray ) * 25.0/self._beta/self._beta

	def verboten_score(self, seq):
		return self.spuriouslib.score_verboten_num( seq.tostring(), self.N, self.wc_aref, self.eq_aref )

	def mutate( self, seq, n ):
		"""Randomly mutate n free bases."""
		
		ind = random.sample(self.freeindex,n)
		
		for i in ind:
			cs = choicearraybin[self.choices[i]]
			new = random.sample(cs,1)[0]
			seq[i] = new
		
		return self.prop( seq )

	def randseq( self ):
		seq = np.zeros(self.N,dtype='int8')

		for i in xrange(0,len(seq)):
			if self.nchoices[i]==0:
				seq[i]=0
				continue
			seq[i] = random.sample(choicearraybin[self.choices[i]],1)[0]
		
		return self.prop(seq)

	def mutate_single(self, seq):
		ind = self.freeindex[random.randint(0,self.lenfree-1)]

		cs = choicearraybin[self.choices[ind]]
		new = random.sample(cs,1)[0]
		seq[ind] = new
		seq[self.constrained & (self.eq == ind+1)] = new
		seq[self.constrained & (self.wc == ind+1)] = 5-new

	def prop( self, seq ):
		"""propagate unconstrained params"""
		seq[self.constrained] = seq[self.eq[self.constrained]-1]
		seq[self.constrained] = 5-seq[self.wc[self.constrained]-1]
		#seq[self.blank] = 0
		return seq

def wcletter(array):
	return np.array([wcletterdict[x] for x in array])

def main():
	import sys
	import logging

	logging.getLogger().setLevel(logging.DEBUG)
	#start = open(sys.argv[1]).read()
	template = open(sys.argv[1]).read()
	wc = np.loadtxt(sys.argv[2])
	eq = np.loadtxt(sys.argv[3])
	if len(sys.argv) == 5:
		time = float(sys.argv[4])
		auto = 1
	else:
		tmax = float(sys.argv[4])
		tmin = float(sys.argv[5])
		steps = float(sys.argv[6])
		updates = float(sys.argv[7])
		auto = 0

	weights = np.array([0.08,0.08,0.08])
	sp_weight = 1.3
	vb_weight = 1.0

	import anneal

	
	system = SpuriousSystemNum(template,wc,eq,spurious_weights=weights)

	sb = system.randseq()
	scorer = lambda x: sp_weight*system.spurious_score(x)+vb_weight*system.verboten_score(x)
	annealer = anneal.Annealer( scorer , system.mutate_single )

	log.info( "Initial sequence:" )
	log.info( ''.join([binbaser[x] for x in sb]) )
	log.info( "" )
	log.info( "Have score of %f. Verboten score is %f. Spurious score is %f." % (scorer(sb),system.verboten_score(sb),system.spurious_score(sb)) )
	log.info( "" )
	log.info( "Spurious Info:" )
	log.info( "" )
	log.info( str(system.spurious(sb)) )
	log.info( "" )
	log.info( str(system.spurious(sb,mismatch=1,kmin=5,kmax=10)) )
	log.info( "" )
	log.info( "Annealing" )

	if auto == 1:
		out = annealer.auto( sb, time )
	else:
		out = annealer.anneal( sb, tmax, tmin, steps, updates )

	log.info( "DONE" )
	log.info( "" )
	log.info( "Have score of %f. Verboten score is %f. Spurious score is %f." % (out[1],system.verboten_score(out[0]),system.spurious_score(out[0])) )
	log.info( "" )
	log.info( "Spurious Info:" )
	log.info( "" )
	log.info( str(system.spurious(out[0])) )
	log.info( "" )
	log.info( str(system.spurious(out[0],mismatch=1,kmin=5,kmax=10)) )
	log.info( "" )
	print ''.join([binbaser[x] for x in out[0]])



	
if __name__ == '__main__':
	main()



		
