/*

  Converted to C from Matlab by hand, EW, 4/00 

  command-line args enhanced, EW, 8/02    

  Merged with spuriousCfold.c, EW, 8/02.  Now *requires* the ViennaRNA package! 

gcc -Wall -O3 spuriousC.c -o spuriousC -I/research/include/ViennaRNA-1.4 -L/research/lib -D VIENNA_DIR="/research/src/ViennaRNA-1.4/" -lRNA -lm

gcc -Wall -O3 spuriousC.c -o spuriousC -I/Users/winfree/DNA_cluster/ViennaRNA-1.4/H -L/Users/winfree/DNA_cluster/ViennaRNA-1.4/lib -D VIENNA_DIR="/Users/winfree/DNA_cluster/ViennaRNA-1.4/" -lRNA -lm

  %%% also, need to update to use ViennaRNA-1.8, which is the latest as of Dec 2008

  %%% as a command line option, should allow to use NUPACK routines... that would be nice...

  %%% or, as a compile-time option, should be able to DISABLE thermodynamic folds, so it can compile stand-alone

*/


#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>


#include <math.h>
#include <ctype.h>
#include "fold.h"
#include "part_func.h"
#include "fold_vars.h"
#include "PS_dot.h"
#include "utils.h"
extern void  read_parameter_file(const char fname[]);

extern double log2(double);

char randbasec(char Stc);
char randbasec_nonidentity(char Sc, char Stc);

#define max(a,b) ((a)>(b)?(a):(b))
#define MAX(a,b) ((a)>(b)?(a):(b))

#define min(a,b) ((a)<(b)?(a):(b))
#define MIN(a,b) ((a)<(b)?(a):(b))
#define templatemismatch printf("Warning: mismatch between template St and sequence S.\n")

// globals:

int debug=0;  int maxN = 100000; int maxPrevented = 10000; int maxTargeted = 10000; int maxSubsequence=10;
int quiet=0, watch=1; int attempt_factor = 4;
double W_verboten=0.0, W_spurious=1.0, W_nofold=0.0, W_struct=0.0, W_prob=0.0, W_bonds=0.0, W_targeted =0; 
double spurious_beta=5.0, spurious_mm=25.0;
double spurious_intraS=25.0/3125, spurious_interS=5.0/3125, spurious_interC=1.0/3125;
double factors_spurious[18];
double targeted_beta=5.0, targeted_mm = 25.0;
double targeted_intraS=25.0/3125, targeted_interS=5.0/3125;
double factors_targeted[12];

int spurious_equality=1, reversible_spurious=1;
char *nofold_include=NULL;
double nofold_mdG=1.0, nofold_mdGG=1.0, nofold_mddG=1.0, nofold_ndG=1.0, nofold_ndGG=1.0, nofold_nddG=1.0;
int N=0,NrS=0,NSt=0,Nwc=0,Neq=0, Npv=0, Ntar=0,  Nonspace=0, Nequiv=0, mmax=0; 
int *testwc, *testeq; 
char *testS, *testSt, *prevented_s;
char *targeted_s1, *targeted_s2, *targeted_pair, *targeted_pair_comp; 
char *targeted_s1t, *targeted_s2t, *targeted_pairt;
char *rS_filename=NULL, *St_filename=NULL, *wc_filename=NULL, *eq_filename=NULL, *S_filename=NULL, *pv_filename=NULL, *tar_filename=NULL;
char **prevented_list;
int *prevented_lengths; 
long int *prevented_violations;
int *prevented_occurrences_S, *prevented_occurrences_St;
int **global_targeted_list;
int **targeted_pair_wc, **targeted_pair_eq;
int *P_rev, *P_hist,  *p_rev;     // persistent resources for spurious_reversible.
char *S1234, *S1234wc;    // persistent resources for spurious_reversible. 

/* 
% ENERGIES OF FORMATION
%
% nearest-neighbor parameters from nnHS.m in DNAdesign
%
% numbers & eqns from SantaLucia, PNAS 1998, 95, 1460--1465
% quoting "unified" params of 
%            Allawi & SantaLucia, Biochemistry 1997, 36, 10581--10594
%
% nnDH(i,j) gives enthalpy for nearest-neighbor interaction i,j.
% e.g.    1,2 =>  5' AC 3'  paired with 3' TG 5'.
%
% nnDH = [ AA, AC, AG, AT;
%          CA, CC, CG, CT;
%          GA, GC, GG, GT;
%          TA, TC, TG, TT ]
*/

//   using A=1, C=2, G=3, T=4,    other=0;

// enthalpy in Kcal/mole       
double nnDH[5][5] = { {   0,    0,    0,    0,    0 },
                      {   0, -7.9, -8.4, -7.8, -7.2 },
                      {   0, -8.5, -8.9,-10.6, -7.8 },
                      {   0, -8.2, -9.8, -8.0, -8.4 },
                      {   0, -7.2, -8.2, -8.5, -7.9 } };

// entropy in cal/mole
double nnDS[5][5] = { {   0,     0,     0,     0,     0 },
                      {   0, -22.2, -22.4, -21.0, -20.4 },
                      {   0, -22.7, -19.9, -27.2, -21.0 },
                      {   0, -22.2, -24.4, -19.9, -22.4 },
                      {   0, -21.3, -22.2, -22.7, -22.2 } };

// free energy
#define nnDG(i,j) (nnDH[i][j] - (temperature+273.15)*nnDS[i][j]/1000.0)

// ---------------------------------------------------------------------------
// writes Sorig into S, as a sequence of 1,2,3,4 for A,C,G,T, or 0 if none-of-the-above
void make1234(char *S, char *Sorig)
{
  int i;  int N; N = strlen(Sorig);
 //% map A,C,G,T to 1,2,3,4
 for (i=0; i<N; i++) { 
   switch (Sorig[i]) {
    case 'A': case 'a': case 1: S[i]=1; break;
    case 'C': case 'c': case 2: S[i]=2; break;
    case 'G': case 'g': case 3: S[i]=3; break;
    case 'T': case 't': case 4: S[i]=4; break;
    default: S[i]=0; 
   }
 }
}

//  ------------------------------------------------------------------------------
void print1234(char *S)
{
  int i, N; N=strlen(S);
  for (i=0; i<N; i++) {
    switch (S[i]) {
     case (1): printf("A"); break;
     case (2): printf("C"); break;
     case (3): printf("G"); break;
     case (4): printf("T"); break;
     case (0): printf(" "); break;
    }
  }
  printf("\n");	
}

// ------------------------------------------------------------------------------
/*
%function C = spurious(S,kmin,kmax, wc,eq)
% C = spurious(S,kmin,kmax, wc,eq)
%
%   S is a sequence of ACTG or 1234 representing the strands in
%        several DNA complexes to be used in an experiment.
%        Strands are separated by ' ' or 0; the strands in separate
%        complexes are separated by '  ' or 0 0.
%
%   eq(i) gives the index of the least base i must be equal to
%         A 0 indicates that no bases must be equal by design 
%   wc(i) gives the index of the least base i must be complementary to
%         A -1 indicates that no base pairing is designed
% 
%   makes a table Cik where 1 <= i <= 3 (intramolec,intracomplex,intercomplex)
%   and kmin <= k <= kmax.  Cik gives the number of undesired WC subsequences
%   of length k which match another subsequence exactly...
%       i=1  intramolecular              (same strand)
%       i=2  intermolecular intracomplex (different strand, same complex)
%       i=3  intermolecular intercomplex (different strand, different complex)
%   Note: counts palindromes, also self-overlapping matches
%
%   This is meant to be a fast fn when compiled, esp. for large S.
%
%   DNAdesign-0.02          Erik Winfree
*/

/* 
 allocates 3 x (kmax-kmin+1) vector [intram 1..n, intrac 1..n, interc 1..n]

 if wc==NULL or eq==NULL, behaves as spurious0.
 if mismatch==1, behaves as spurious1.  i.e. counts almost-matches containing a single mismatch, ignores wc & eq.
 if spurious_equality==1, then it also counts identity matches, not just WC matches. (only for mismatch==0)
*/

int *spurious(int mismatch, char *Sorig, int kmin, int kmax, int *wc, int *eq)
{
 int *C; char *S,*Swc;
 int i, imax = (kmax-kmin+1); // note pwkr 02/2010, imax also used to be used differently in main, now smax
 int N; int k,j;

 C = calloc(3*imax, sizeof(int)); if (C==NULL) return NULL;
 for (i=0; i<3*imax; i++) C[i] = 0;

 N = strlen(Sorig);  // note: wc & eq must be length N, or ELSE!!! 

 S   = malloc(N); if (S  ==NULL) { free(C); return NULL; }
 Swc = malloc(N); if (Swc==NULL) { free(C); free(S); return NULL; }

 make1234(S,Sorig);    // convert ACGT into 1234

 for (i=0; i<N; i++) {                   //% WC bases, same order
   if (S[i]==0) Swc[i]=0; else Swc[i]=5-S[i]; 
 }

 for (k=kmin; k<=kmax; k++) { 
  int M, *P, *p, Cs, Ss, m, mWC;
  M = 1<<(2*k);  //% 4^k;
  P = calloc(M,sizeof(int));  //% hold index of most-recently-found instance of subsequence m
  p = calloc(N,sizeof(int));  //% holds pointers to previous instance of subsequence at i
  Cs = 0;             //% the start of most recent complex
  Ss = 0;             //% the start of most recent strand
  m = 0;              //% the ID of current subsequence, binary coding of current subsequence
  mWC = 0;            //% the ID for this subsequence's WC complement
  for (i=1; i<=N; i++) {     //  i indexes character position in the sequence string, starting at 1
   if ( S[i-1] == 0 ) {      // detect space,  reset strand & complex indicators
    Ss = i; m=0; mWC = 0;
    if ( i>1 && S[i-1 -1] == 0 ) Cs = i;  // two spaces, complex
   } else {                  // within strand/complex
    m   = ( 4*m + S[i-1]-1 ) % M;    // current subsequence
    mWC = (mWC/4) + (Swc[i-1]-1)*M/4; 
    if (i-Ss >= k && mismatch==0) {  // begin spurious0
       p[i-1]=P[m]; P[m]=i; //% add site to list of subsequences
       j=P[mWC];     //% look for previous WC occurrences, including self
       while (j != 0) { int msum,jjj;
         if (wc != NULL && eq != NULL) { //% exclude expected WC match
           for (msum = 0, jjj=0; jjj<k; jjj++)
              if (wc[i-k+1 +jjj-1] == eq[j -jjj-1]) msum++;
         } else msum=0;
         if (msum != k) {                      
           if (j>Ss) C[       k-kmin+1 -1]++;
           if (j>Cs) C[imax  +k-kmin+1 -1]++;
	             C[2*imax+k-kmin+1 -1]++;
         }
         j = p[j -1];
       } // previous WC occurrences
       if (spurious_equality) {
        j=p[i-1];     //% look for previous equality occurrences, not including self
        while (j != 0) { int msum,jjj;
         if (wc != NULL && eq != NULL) { //% exclude expected EQ match
           for (msum = 0, jjj=0; jjj<k; jjj++)
              if (eq[i -jjj-1] == eq[j -jjj-1]) msum++;
         } else msum=0;
         if (msum != k) {                      
           if (j>Ss) C[       k-kmin+1 -1]++;
           if (j>Cs) C[imax  +k-kmin+1 -1]++;
	             C[2*imax+k-kmin+1 -1]++;
         }
         j = p[j -1];
        }
       } // previous EQ occurrences
    }    // end spurious0
    if (i-Ss >= k && mismatch==1) { int nt, bp, misWC;   // begin spurious1
     for (nt=2; nt<=k-1; nt++) {  //% where is mismatch
      //% misWC = mWC - 4^(nt-1) * rem(fix(mWC/4^(nt-1)),4);
      misWC = mWC & (-1 - (3<<(2*nt-2)));
      for (bp=0; bp<4; bp++) { //% what is value at mismatch
       misWC += bp<<(2*nt-2);
       j=P[misWC];   //% look for previous occurrences, including self
       while (j != 0 && misWC != mWC) { //% count all mismatches only
         if (j>Ss) C[       k-kmin+1 -1]++;
         if (j>Cs) C[imax  +k-kmin+1 -1]++;
            C[2*imax+k-kmin+1 -1]++;
         j = p[j -1];
       }
       if (m==misWC && m!=mWC) { //% "palindrome-with-mismatch"
          C[       k-kmin+1 -1]++;
          C[imax  +k-kmin+1 -1]++;
          C[2*imax+k-kmin+1 -1]++;
       }
       misWC -= bp<<(2*nt-2);
      }
     }
     p[i -1]=P[m]; P[m]=i; //% add this site to list of subsequences
    }  // end spurious 1
   }   // within complex
  }    // for all i positions in S
  C[2*imax+k-kmin+1 -1] = C[2*imax+k-kmin+1 -1] - C[imax+k-kmin+1 -1];
  C[  imax+k-kmin+1 -1] = C[  imax+k-kmin+1 -1] - C[     k-kmin+1 -1];
  free(p); free(P);
 }
 free(S); free(Swc);
 return C;
}


// ------------------------------------------------------------------------------
// This function is the same as spurious but it uses globally defined arrays for the hash table
// and other resources, such as a numerical encoding of the sequence
// rather than re-allocating them every time it is called. 
// It cleans up the hash-table in time linear on the length of the input sequence.
//
// S is over 1234, Swc is the complement over 1234.
int *spurious_reversible(int mismatch, char *S, char *Swc, int kmin, int kmax, int *wc, int *eq)
{
 int *C; 
 int i,ii, imax = (kmax-kmin+1); // note pwkr 02/2010, imax also used to be used differently in main, now smax
 int N; int k,j, pcount=0;

 C = calloc(3*imax, sizeof(int)); if (C==NULL) return NULL;
 for (i=0; i<3*imax; i++) C[i] = 0;

 N = strlen(S);  // note: wc & eq must be length N, or ELSE!!! 

 for (k=kmin; k<=kmax; k++) { 
  int M, Cs, Ss, m, mWC;
  M = 1<<(2*k);  //% 4^k;
  // P_rev and p_rev now global. 
  // P_rev = calloc(M,sizeof(int));  //% hold index of most-recently-found instance of subsequence m
  // p_rev = calloc(N,sizeof(int));  //% holds pointers to previous instance of subsequence at i
  Cs = 0;             //% the start of most recent complex
  Ss = 0;             //% the start of most recent strand
  m = 0;              //% the ID of current subsequence, binary coding of current subsequence
  mWC = 0;            //% the ID for this subsequence's WC complement
  pcount = 0;         //P_rev and P_hist are clean.
  for (i=1; i<=N; i++) {     //  i indexes character position in the sequence string, starting at 1
   if ( S[i-1] == 0 ) {      // detect space,  reset strand & complex indicators
    Ss = i; m=0; mWC = 0;
    if ( i>1 && S[i-1 -1] == 0 ) Cs = i;  // two spaces, complex
   } else {                  // within strand/complex
    m   = ( 4*m + S[i-1]-1 ) % M;    // current subsequence
    mWC = (mWC/4) + (Swc[i-1]-1)*M/4; 
    if (i-Ss >= k && mismatch==0) {  // begin spurious0
      p_rev[i-1]=P_rev[m]; P_rev[m]=i; P_hist[pcount] = m; pcount++;//% add site to list of subsequences
       j=P_rev[mWC];     //% look for previous WC occurrences, including self
       while (j != 0) { int msum,jjj;
         if (wc != NULL && eq != NULL) { //% exclude expected WC match
           for (msum = 0, jjj=0; jjj<k; jjj++)
              if (wc[i-k+1 +jjj-1] == eq[j -jjj-1]) msum++;
         } else msum=0;
         if (msum != k) {                      
           if (j>Ss) C[       k-kmin+1 -1]++;
           if (j>Cs) C[imax  +k-kmin+1 -1]++;
	             C[2*imax+k-kmin+1 -1]++;
         }
         j = p_rev[j -1];
       } // previous WC occurrences
       if (spurious_equality) {
        j=p_rev[i-1];     //% look for previous equality occurrences, not including self
        while (j != 0) { int msum,jjj;
         if (wc != NULL && eq != NULL) { //% exclude expected EQ match
           for (msum = 0, jjj=0; jjj<k; jjj++)
              if (eq[i -jjj-1] == eq[j -jjj-1]) msum++;
         } else msum=0;
         if (msum != k) {                      
           if (j>Ss) C[       k-kmin+1 -1]++;
           if (j>Cs) C[imax  +k-kmin+1 -1]++;
	             C[2*imax+k-kmin+1 -1]++;
         }
         j = p_rev[j -1];
        }
       } // previous EQ occurrences
    }    // end spurious0
    if (i-Ss >= k && mismatch==1) {  // begin spurious1
      int nt, bp, misWC;   
      for (nt=2; nt<=k-1; nt++) {  //% where is mismatch
	//% misWC = mWC - 4^(nt-1) * rem(fix(mWC/4^(nt-1)),4);
	misWC = mWC & (-1 - (3<<(2*nt-2)));
	for (bp=0; bp<4; bp++) { //% what is value at mismatch
	  misWC += bp<<(2*nt-2);
	  j=P_rev[misWC];   //% look for previous occurrences, including self
	  while (j != 0 && misWC != mWC) { //% count all mismatches only
	    if (j>Ss) C[       k-kmin+1 -1]++;
	    if (j>Cs) C[imax  +k-kmin+1 -1]++;
            C[2*imax+k-kmin+1 -1]++;
	    j = p_rev[j -1];
	  }
	  if (m==misWC && m!=mWC) { //% "palindrome-with-mismatch"
	    C[       k-kmin+1 -1]++;
	    C[imax  +k-kmin+1 -1]++;
	    C[2*imax+k-kmin+1 -1]++;
	  }
	  misWC -= bp<<(2*nt-2);
	}
      }
      p_rev[i -1]=P_rev[m]; P_rev[m]=i; P_hist[pcount] = m; pcount++;//% add this site to list of subsequences
    }  // end spurious 1
   }   // within complex
  }    // for all i positions in S
  C[2*imax+k-kmin+1 -1] = C[2*imax+k-kmin+1 -1] - C[imax+k-kmin+1 -1];
  C[  imax+k-kmin+1 -1] = C[  imax+k-kmin+1 -1] - C[     k-kmin+1 -1];
  // cleanup p_rev and P_rev
  for (ii=0; ii< N; ii++)  p_rev[ii] = 0;
  for (ii=0; ii<pcount; ii++) P_rev[P_hist[ii]] = 0;  
 } // for all values of subsequence size k.
 return C;

}




// -------------------------------------------------------------------------------

// c MUST BE UPPERCASE ACGT & degenerate bases
char WC(char c)
{
    switch(c) {
    case 'A': return 'T';
    case 'C': return 'G';
    case 'G': return 'C';
    case 'T': return 'A';
    case 'R': return 'Y';
    case 'Y': return 'R';
    case 'W': return 'W';
    case 'S': return 'S';
    case 'K': return 'M';
    case 'M': return 'K';
    case 'B': return 'V';
    case 'V': return 'B';
    case 'D': return 'H';
    case 'H': return 'D';
    case 'N': return 'N';
    default: return ' '; 
    }
}

// ---------------------------------------------------------------------------------
// modifies S.   S MUST BE UPPERCASE ACGT & degenerate bases
// (Note: pwkr, 02/2010 this function is never used, and doesn't create a reverse complement)
void WCstring(char *S)
{
  int i,n; n=strlen(S);
  for (i=0; i<n; i++) { S[i]=WC(S[i]); }
}

// ------------------------------------------------------------------------------------
// Takes a sequence over 1234 and makes the WC complement, but doesn't change the order
// so the result Swc is 3' to 5'
void WCstring1234(char *S, char *Swc)
{
  int i, N;
  N = strlen(S);
  for (i=0; i<N; i++) {                   //% WC bases, same order
    if (S[i]==0) Swc[i]=0; else Swc[i]=5-S[i]; 
  }
}

// ---------------------------------------------------------------------------------
// chars c1 must be uppercase ACGT & degenerate bases, and c2 must be uppercase ACGT
// this function returns true if c2 is a base, A,C,G, or T that is compatible with c1
// note that this function is not commutative. 
// follows strcmp convention, compatibility returns 0
int compatible_base_cmp(char c1, char c2)
{
  switch (c1) {
  case 'A': switch (c2) {case 'A': return 0; default: return 1;}
  case 'C': switch (c2) {case 'C': return 0; default: return 1;}
  case 'G': switch (c2) {case 'G': return 0; default: return 1;}
  case 'T': switch (c2) {case 'T': return 0; default: return 1;}
  case 'R': switch (c2) {case 'A': case 'G': return 0; default: return 1;}
  case 'Y': switch (c2) {case 'C': case 'T': return 0; default: return 1;}
  case 'W': switch (c2) {case 'A': case 'T': return 0; default: return 1;}
  case 'S': switch (c2) {case 'C': case 'G': return 0; default: return 1;}
  case 'K': switch (c2) {case 'G': case 'T': return 0; default: return 1;}
  case 'M': switch (c2) {case 'A': case 'C': return 0; default: return 1;}
  case 'B': switch (c2) {case 'C': case 'G': case 'T': return 0; default: return 1;}
  case 'V': switch (c2) {case 'A': case 'C': case 'G': return 0; default: return 1;}
  case 'D': switch (c2) {case 'A': case 'G': case 'T': return 0; default: return 1;}
  case 'H': switch (c2) {case 'A': case 'C': case 'T': return 0; default: return 1;}
  case 'N': switch (c2) {case 'A': case 'C': case 'G': case 'T': return 0; default: return 1;}
  case ' ': switch (c2) {case ' ': return 0; default: return 1;}
  default: return -1; // error condition, also will happen when there are spaces...
  }
}

// ---------------------------------------------------------------------------------
// strings s1 must be uppercase ACGT & degenerate bases, and s2 must be uppercase ACGT
// this function returns true if s2 is a base, A,C,G, or T that is compatible with s1
// note that this function is not commutative. 
// follows strcmp convention, compatibility returns 0
int compatible_sequence_cmp(char *s1, char *s2)
{
  int i, slen1, slen2, num_compatible = 0;
  slen1 = strlen(s1);
  slen2 = strlen(s2);
  if (slen1 != slen2) return -1; // strings unequal length
  for (i=0; i < slen1; i++) {
    if (!compatible_base_cmp(s1[i],s2[i])) num_compatible++;
  }
  if (num_compatible == slen1) return 0; // sequences are compatible
  else return 1; // sequences are the same length but not compatible.
}

// ---------------------------------------------------------------------------------
int count_nonspace(char *St)
{
  int i, num_space=0, slen;
  slen = strlen(St);
  for (i=0; i<slen; i++) {
    switch (St[i]) {
    case ' ': num_space++;
    default: ;
    }
  }
  return (slen - num_space);
}
// ---------------------------------------------------------------------------------
// counts the number of bases fixed as A, G, C, and T in sequence template. 
int count_fixed(char *St)
{
  int i, num_fixed=0, slen;
  slen = strlen(St);
  for (i=0; i<slen; i++) {
    switch (St[i]) {
    case 'A': case 'G': case 'C': case 'T': num_fixed++;
    default: ;
    }
  }
  return num_fixed;
}

// ---------------------------------------------------------------------------------
// returns zero if two arrays indicating the number of prevented sequences of each type
// for two different strings (ie. testS and testSt) are the same. Returns the sum
// of their differences otherwise. 
int prevent_cmp(int *pv1, int *pv2, int num_pv)
{
  int i, prevent_differences=0;
  for (i=0; i<num_pv; i++) {
    prevent_differences += abs(pv1[i] - pv2[i]);
  }
  return prevent_differences;
}

// ---------------------------------------------------------------------------------
// takes a sequence S over ACGT an a template St over ACGT+degenerate
// and returns a count of the places S is inconsistent with St. 
// Modifies S to be consistent with St.
int template_match(char *S, char *St)
{
  int n, i, num_matches=0;
  n = strlen(S);
  for (i=0; i< n; i++) {
    if (compatible_base_cmp(St[i],S[i])==0) // seed is compatible with template at i
      num_matches++;
    else {
      S[i] = randbasec(St[i]);   // force seed to be compatible with template at i
    }
  }
  return num_matches;
}

// --------------------------------------------------------------------------
/*
% function S = constrain(rS, wc, eq)
% S = constrain(rS, wc, eq)
%
%   rS is the raw sequence of 'ACGT '; 
%      bases that should be WC or equal might not be. 
%   wc and eq are as built by constraints().
%
%   S is rS, but with each equivalence class's rep's bases used to
%            determine the entire class and its WC class,
%            so now wc and eq are obeyed.
%
%   DNAdesign-0.02          Erik Winfree
*/

/*
   modifies S.  S MUST BE UPPERCASE ACGT & degenerate bases
*/
void constrain(char *S, int *wc, int *eq)
{
  int n,i,j; int *marked;

  n=strlen(S);
  marked = (int *) calloc(n, sizeof(int));

  for (i=0; i<n; i++) {
    if (!marked[i]) {
      // class = find(eq==eq(i));
      // for j=1:length(class), S(class(j))=S(i); end
      // marked(class)=1;

      for (j=0; j<n; j++) {
        if (eq[j] == eq[i]) { S[j] = S[i]; marked[j]=1; }
      }

      // class = find(eq==wc(i));
      // for j=1:length(class), S(class(j))=WC(S(i)); end
      // marked(class)=1;

      for (j=0; j<n; j++) {
        if (eq[j] == wc[i]) { S[j] = WC(S[i]); marked[j]=1; }
      }
    }
  }
  free(marked);

}
// --------------------------------------------------------------------------

/*
% function bad6mer = make_bad6mer()
% bad6mer = make_bad6mer()
% 
%   Make a table of costs for "verboten" 6-base subsequences.
%   The costs are just made up... vaguely based on rules of thumb
%   such as "don't have 3 Gs in a row"... "avoid many purines in
%   a row; you might get a triple helix"..."don't have long AT runs"...
%   etc...   Some may by myths or mistakes.
%
%   The results are saved in the file Bad6mer.mat,
%   so it only needs to be run once.
%
%   The only reason this is a function is so it can be compiled, e.g.
%      mcc -wirh make_bad6mer
%
%   DNAdesign-0.02          Erik Winfree
*/

double bad6mer[1<<12];

void make_bad6mer()   
{
  char S[7]; char ACGT[5]=" ACGT";
  int i,j, b1,b2,b3,b4,b5,b6;
  int co[5],ce[5],c[5];
  double ATpenalty = 0.1;

  for (i=0; i < (1<<12);  i++) bad6mer[i] = 0.0;

  i=0;
  for (b1=1; b1<=4; b1++)
   for (b2=1; b2<=4; b2++)
    for (b3=1; b3<=4; b3++)
     for (b4=1; b4<=4; b4++)
      for (b5=1; b5<=4; b5++)
       for (b6=1; b6<=4; b6++) {
         S[0]=ACGT[b1]; S[1]=ACGT[b2]; S[2]=ACGT[b3];
         S[3]=ACGT[b4]; S[4]=ACGT[b5]; S[5]=ACGT[b6]; S[6]=0;
         //% c(i) = how many base type i are present.  
         //% count c(i) for even & odd positions
         for (j=1;j<=4;j++) {co[j]=ce[j]=0;}
         co[b1]=co[b1]+1; ce[b2]=ce[b2]+1; co[b3]=co[b3]+1;
         ce[b4]=ce[b4]+1; co[b5]=co[b5]+1; ce[b6]=ce[b6]+1;
         for (j=1;j<=4;j++) {c[j]=co[j]+ce[j];}
         bad6mer[i] = 
	  (c[2]+c[3]==0) +             //% all AT
	  (c[1]+c[4]==0) +             //% all CG
	  (c[1]+c[3]==0) +             //% all AG (Pu)
	  (c[2]+c[4]==0) +             //% all TC (Py)
	  (co[1]+co[3]+ce[2]+ce[4]==0) +  //% alt PuPy
	  (ce[1]+ce[3]+co[2]+co[4]==0) +  //% alt PyPu
          (strstr(S,"GGG")!=NULL) + 
          (strstr(S,"GGGG")!=NULL) + 
          (strstr(S,"CCC")!=NULL) + 
          (strstr(S,"CCCC")!=NULL) + 
          (strstr(S,"TTTT")!=NULL) + 
          (strstr(S,"AAAA")!=NULL) + 
          (strstr(S,"GCGC")!=NULL) + 
          (strstr(S,"GGCC")!=NULL) + 
          (strstr(S,"CCGG")!=NULL) + 
          (strstr(S,"CGCG")!=NULL) + 
	  ((c[1]+c[4])==5 && (b1==2 || b1==3)) +  //% 5 AT is bad
	  ((c[1]+c[4])==5 && (b6==2 || b6==3)) +  //% 5 AT is bad
	  ((c[2]+c[3])==5 && (b1==1 || b1==4)) +  //% 5 CG is bad
	  ((c[2]+c[3])==5 && (b6==1 || b6==4)) +  //% 5 CG is bad
	  ATpenalty*(c[1]+c[4]>3);            //% try to increase %CG
         i=i+1;
       }
}


// --------------------------------------------------------------------------
/*
%function V = verboten(S)
% V = verboten(S)
%   counts number of occurrences of to-be-avoided sequences, e.g.
%      6-mers containing GGG GGGG TTTT AAAA CCCC GCGC GGCC CCGG CGCG 
%      6-mer tracts of all-AT, all-GC, all-Pu, all-Py, and alternating Pu Py
%
%   NOTE: the variable 'bad6mer' must be defined as global to indicate
%   penalties; it should loaded from Bad6mer.mat, in turn created once
%   by make_bad6mer().
%
%   DNAdesign-0.02          Erik Winfree
*/

double score_verboten(char *Sorig, int *wc, int *eq)
{
 int i, K, N, M, m; char *S; double V;

 N = strlen(Sorig);  
 S = malloc(N); if (S==NULL) { return 0.0; }

 //% map A,C,G,T to 1,2,3,4
 for (i=0; i<N; i++) { 
   switch (Sorig[i]) {
    case 'A': case 'a': case 1: S[i]=1; break;
    case 'C': case 'c': case 2: S[i]=2; break;
    case 'G': case 'g': case 3: S[i]=3; break;
    case 'T': case 't': case 4: S[i]=4; break;
    default: S[i]=0; 
   }
 }
 V=0.0; K=5;

 M = 1<<12;         // M = 10000000000000
 m = 0;                //% the ID of current subsequence
 for (i=0; i<N; i++) {
  if (S[i] == 0) {     //% reset strand & complex indicators
    m=0; K=5;          //% wait till full 6-mer considered
  } else {
    m = (  (m<<2) + S[i]-1 ) & (M-1);   //shift 2 bit char left, add new 2 bits of current char, mask off
                                        //12 bits with 111111111111 to kill left two bits.
    if (K) K=K-1; else V=V+bad6mer[m]; 
  }				       
 }


// printf("verboten('%s')=%f\n\n",Sorig,V);

 free(S);
 return V; 
}

//---------------------------------------------------------------------------------
int *targeted_spurious(int mismatch, char *S, int kmin, int kmax,
		       int *wc, int *eq, int num_tar, int **targeted_list)
{
  int i,j, s1s, s1e, s2s, s2e, s1len, s2len, spairlen; 
  // s1s is string 1 start, s1e is string 2 end, etc for s2s, s2e
  int *C,*C_sum;
  int imax = (kmax-kmin+1);

  C_sum = calloc(3*imax, sizeof(int)); if (C_sum==NULL) return NULL;
  for (i=0; i<3*imax; i++) C_sum[i] = 0;

  for (i=0; i<num_tar; i++) {
    // get indices of pair i of targeted strings,
    // they start with 1 so  subtract 1 so they match C array convention.
    s1s = targeted_list[i][0]-1;
    s1e = targeted_list[i][1]-1;
    s2s = targeted_list[i][2]-1;
    s2e = targeted_list[i][3]-1;
    s1len = s1e - s1s + 1;
    s2len = s2e - s2s + 1;
    // copy targeted strings out of S, to targeted_s1, targeted_s2, global strings of size maxTargeted
    strncpy(targeted_s1, S + s1s, s1len);
    strncpy(targeted_s2, S + s2s, s2len);
    // null terminate
    targeted_s1[s1len] = '\0';
    targeted_s2[s2len] = '\0';
    // catenate to create a two strand system
    strcpy(targeted_pair, targeted_s1);
    strcat(targeted_pair, " ");
    strcat(targeted_pair, targeted_s2);
    spairlen = strlen(targeted_pair);
    if (spairlen != s1len + s2len + 1) printf("Warning: error catenating targets for targeted spurious.\n");

    // call spurious accumulate in C_sum
    C = spurious(mismatch, targeted_pair, kmin, kmax, targeted_pair_wc[i], targeted_pair_eq[i]);
    for (j=0; j<3*imax; j++) C_sum[j] += C[j];
  } // for i over all elements of the target list. 
  return C_sum;
}

//---------------------------------------------------------------------------------
//  S and Swc are over 1234 but they are still strings...
int *targeted_spurious_reversible(int mismatch, char *S, char *Swc, int kmin, int kmax,
				  int *wc, int *eq, int num_tar, int **targeted_list)
{
  int i,j, s1s, s1e, s2s, s2e, s1len, s2len, spairlen; 
  // s1s is string 1 start, s1e is string 2 end, etc for s2s, s2e
  int *C,*C_sum;
  int imax = (kmax-kmin+1);

  C_sum = calloc(3*imax, sizeof(int)); if (C_sum==NULL) return NULL;
  for (i=0; i<3*imax; i++) C_sum[i] = 0;

  for (i=0; i< num_tar; i++){
    // get indices of pair i of targeted strings,
    // they start with 1 (from matlab) so  subtract 1 so they match C array convention.
    s1s = targeted_list[i][0]-1;
    s1e = targeted_list[i][1]-1;
    s2s = targeted_list[i][2]-1;
    s2e = targeted_list[i][3]-1;
    s1len = s1e - s1s + 1;
    s2len = s2e - s2s + 1;
    // copy targeted strings out of S, to targeted_s1, targeted_s2, global strings of size maxTargeted
    strncpy(targeted_s1, S + s1s, s1len);
    strncpy(targeted_s2, S + s2s, s2len);
    // null terminate
    targeted_s1[s1len] = '\0';
    targeted_s2[s2len] = '\0';
    // catenate to create a two strand system
    strcpy(targeted_pair, targeted_s1);
    strcat(targeted_pair, " ");
    strcat(targeted_pair, targeted_s2);
    spairlen = strlen(targeted_pair);
    if (spairlen != s1len + s2len + 1) printf("Warning: error catenating targets for targeted spurious.\n");
    // make targeted_pair_comp the watson-crick of targeted_pair
    WCstring1234(targeted_pair, targeted_pair_comp);
    targeted_pair_comp[strlen(targeted_pair)] = '\0'; // null terminate
    if (strlen(targeted_pair) != strlen(targeted_pair_comp)) {
      printf("Warning, screwed up null termination for target system.\n");
    }
    else {
      printf("all ok with string termination, remove check.\n");
    }

    // call spurious accumulate in C_sum
    C = spurious_reversible(mismatch, targeted_pair, targeted_pair_comp,
			    kmin, kmax, targeted_pair_wc[i], targeted_pair_eq[i]);
    for (j=0; j<3*imax; j++) C_sum[j] += C[j];
  } // for i over the target list
  return C_sum;
}
//---------------------------------------------------------------------------------
// function score_targeted(char *S, int *wc, int *eq, int num_tar, int **targeted_list)
//
// This function takes a list of pairs of substrings of testS, and for
// each pair performs a spurious test on that pair. 
// It adds these spurious tests up and returns the score. 
// The point is that while the normal spurious score takes into
// account spurious mismatches for all possible subsequences of S, 
// it may be that certain subsequences of S are more prone to bind to 
// eachother than others, because, for example, they are held in close
// proximity by specific structure. If we want these sequences to _not_ bind
// we may wish to increase the relative weight of their spurious score. (pwkr, 02/2010)
// it is OK for a pair of subsequences to have wc or eq matches to eachother _but_
// wc or eq matches to other substrings will be lost. 
// assumes substring 1 indices < substring 2 (enforced on loading)

double score_targeted(char *S, int *wc, int *eq, int num_tar, int **targeted_list)
{
  int i;
  int *C; 
  double score=0.0; 
 
  C = targeted_spurious(0,S,3,8,wc,eq,num_tar,targeted_list);
  for (i=0; i<12; i++) score+=factors_targeted[i]*C[i];

  C = targeted_spurious(1,S,5,10,NULL,NULL,num_tar,targeted_list);
  for (i=0; i<12; i++) score+=factors_targeted[i]*C[i]*targeted_mm/targeted_beta/targeted_beta;

  return score;
}


// ------------------------------------------------------------------
// score_targeted for use with spurious_reversible
// S and Swc over 1234
double score_targeted_reversible(char *S, char *Swc, int *wc, int *eq, int num_tar, int **targeted_list)
{
  int i;
  int *C; 
  double score=0.0; 
 
  C = targeted_spurious_reversible(0,S,Swc,3,8,wc,eq,num_tar,targeted_list);
  for (i=0; i<12; i++) score+=factors_targeted[i]*C[i];

  C = targeted_spurious_reversible(1,S,Swc,5,10,NULL,NULL,num_tar,targeted_list);
  for (i=0; i<12; i++) score+=factors_targeted[i]*C[i]*targeted_mm/targeted_beta/targeted_beta;

  return score;

}


//---------------------------------------------------------------------------------
// optimizing score_spurious performs so poorly on Milo's test;
// Niles & Robert's positive+negative approach is more successful...
// add some element of positive design...
//    here, we sum the nearest-neighbor Hbond/stacking energy for
//    all desired WC interactions, in Kcal/mole. 
double score_bonds(char *Sorig, int *wc, int *eq)
{
  double score=0; int i,N; char *S;

  //  for (i=0; i<5; i++) 
  //     { for (N=0; N<5; N++) printf("%4.2f ",nnDH[i][N]); printf("\n"); }

  // temperature = 37; 

  N=strlen(Sorig); S = malloc(N); make1234(S,Sorig);

  for (i=1; i<N; i++)
    if (wc[i-1]==wc[i]+1) score += nnDG(S[i-1]+0,S[i]+0);

  return score/2;  // every nt involved in WC get counted twice per stack 
}

//---------------------------------------------------------------
void make_upper(char *S)
{ int i;
  for (i=0; i< strlen(S); i++)  S[i] = toupper(S[i]);
}

// --------------------------------------------------------------------------
/*
%                 1 2 3 4 5  6  7  8  9  10 11  12  13  14  15   0=space
%   degenerates:  A C G T R  Y  W  S  M  K  B   D   H   V   N 
%                         AG CT AT CG AC GT CGT AGT ACT ACG ACGT 
*/
char randbasec(char Stc) 
{
  char *choices;
  switch(toupper(Stc)) { 
  case 'A': choices = "A"; break;
  case 'C': choices = "C"; break;
  case 'G': choices = "G"; break;
  case 'T': choices = "T"; break;
  case 'R': choices = "AG"; break;
  case 'Y': choices = "CT"; break;
  case 'W': choices = "AT"; break;
  case 'S': choices = "CG"; break;
  case 'M': choices = "AC"; break;
  case 'K': choices = "GT"; break;
  case 'B': choices = "CGT"; break;
  case 'D': choices = "AGT"; break;
  case 'H': choices = "ACT"; break;
  case 'V': choices = "ACG"; break;
  case 'N': choices = "ACGT"; break;
  case ' ': choices = " "; break;
  default: printf("Warning: invalid letter, non-space base in template St.\n"); choices = " ";
  }
  return choices[ int_urn(0,strlen(choices)-1) ];
}

// --------------------------------------------------------------------------

/*
%                 1 2 3 4 5  6  7  8  9  10 11  12  13  14  15   0=space
%   degenerates:  A C G T R  Y  W  S  M  K  B   D   H   V   N 
%                         AG CT AT CG AC GT CGT AGT ACT ACG ACGT 
*/
// choses a random base of type Stc that is not equal to Sc
char randbasec_nonidentity(char Sc, char Stc) 
{
  char *choices;

  switch (Sc) {
  case 'A':
    switch (Stc) { 
    case 'A': choices = "A"; break;
    case 'C': templatemismatch; choices = "C"; break;
    case 'G': templatemismatch; choices = "G"; break;
    case 'T': templatemismatch; choices = "T"; break;
    case 'R': choices = "G"; break;
    case 'Y': templatemismatch; choices = "CT"; break;
    case 'W': choices = "T"; break;
    case 'S': templatemismatch; choices = "CG"; break;
    case 'M': choices = "C"; break;
    case 'K': templatemismatch; choices = "GT"; break;
    case 'B': templatemismatch; choices = "CGT"; break;
    case 'D': choices = "GT"; break;
    case 'H': choices = "CT"; break;
    case 'V': choices = "CG"; break;
    case 'N': choices = "CGT"; break;
    default: templatemismatch;  choices = " ";
    } break;
  case 'C':
    switch(Stc) { 
    case 'A': templatemismatch; choices = "A"; break;
    case 'C': choices = "C"; break;
    case 'G': templatemismatch; choices = "G"; break;
    case 'T': templatemismatch; choices = "T"; break;
    case 'R': templatemismatch; choices = "AG"; break;
    case 'Y': choices = "T"; break;
    case 'W': templatemismatch; choices = "AT"; break;
    case 'S': choices = "G"; break;
    case 'M': choices = "A"; break;
    case 'K': templatemismatch; choices = "GT"; break;
    case 'B': choices = "GT"; break;
    case 'D': templatemismatch; choices = "AGT"; break;
    case 'H': choices = "AT"; break;
    case 'V': choices = "AG"; break;
    case 'N': choices = "AGT"; break;
    default: templatemismatch;  choices = " ";
    } break;
  case 'G':
    switch(Stc) { 
    case 'A': templatemismatch; choices = "A"; break;
    case 'C': templatemismatch; choices = "C"; break;
    case 'G': choices = "G"; break;
    case 'T': templatemismatch; choices = "T"; break;
    case 'R': choices = "A"; break;
    case 'Y': templatemismatch; choices = "CT"; break;
    case 'W': templatemismatch; choices = "AT"; break;
    case 'S': choices = "C"; break;
    case 'M': templatemismatch; choices = "AC"; break;
    case 'K': choices = "T"; break;
    case 'B': choices = "CT"; break;
    case 'D': choices = "AT"; break;
    case 'H': templatemismatch; choices = "ACT"; break;
    case 'V': choices = "AC"; break;
    case 'N': choices = "ACT"; break;
    default:  templatemismatch; choices = " ";
    } break;
  case 'T':
    switch(Stc) { 
    case 'A': templatemismatch; choices = "A"; break;
    case 'C': templatemismatch; choices = "C"; break;
    case 'G': templatemismatch; choices = "G"; break;
    case 'T': choices = "T"; break;
    case 'R': templatemismatch; choices = "AG"; break;
    case 'Y': choices = "C"; break;
    case 'W': choices = "A"; break;
    case 'S': templatemismatch; choices = "CG"; break;
    case 'M': templatemismatch; choices = "AC"; break;
    case 'K': choices = "G"; break;
    case 'B': choices = "CG"; break;
    case 'D': choices = "AG"; break;
    case 'H': choices = "AC"; break;
    case 'V': templatemismatch; choices = "ACG"; break;
    case 'N': choices = "ACG"; break;
    default:  templatemismatch; choices = " "; 
    } break;
  default: choices = " ";
  }
  return choices[ int_urn(0,strlen(choices)-1) ];
}

//------------------------------------------------------------------------
// make a single random change, apply the constraints, and insist that the change wasn't eliminated
int mutate(char *S, char *St, int *wc, int *eq, char **pv_list, int *pv_lengths, int num_pv, long int *pv_violations) 
{
  int i,j,k, slen, pvlen, attempts=0; char oldc; 
  slen = strlen(St);
  do {
    attempts++;
    do {
      i = int_urn(0,slen-1);  // generate random base position (nonspace)
    } while (S[i] == ' ');
    oldc = S[i];
    S[i] = randbasec_nonidentity(S[i], St[i]); // base not equal to S[i]
    constrain(S,wc,eq);     // may revert base if it violates wc or eq  constraints
    // now check to see if a prevented sequence occurs in a place
    // where it shouldn't.
    if (oldc != S[i]) {
      for (j=0; j < num_pv; j++) {  // iterate over prevented strings
	pvlen = strlen(pv_list[j]);
	for (k=max(0,i-(pvlen-1)); k<min(i+1,slen-(pvlen-1)); k++) {
	  strncpy(prevented_s, S + k, pvlen);  // prevented_s is a global, maxPrevented size char array
	  prevented_s[pvlen] = '\0'; // null terminate prevented_s
	  if (debug) printf("pvlen %d k %d prevent list %s, %d matched to substring %s, %d, strcmp %d, comp_seq %d.\n", pvlen, k,
			    pv_list[j], strlen(pv_list[j]), prevented_s, strlen(prevented_s), 
			    strcmp(pv_list[j],prevented_s), compatible_sequence_cmp(pv_list[j],prevented_s));
	  if (compatible_sequence_cmp(pv_list[j],prevented_s)==0) // prevented string detected
	    {S[i] = oldc; pv_violations[j]++;}         // Restore S, setup to randomize again
	}
      }
    }
  } while ((oldc == S[i]) && (attempts < mmax));
  // if no base change allowed by template could alter constrained S, this might never halts
  // hence the extra stop condition with 'attempts'.
  return attempts;
}

// --------------------------------------------------------------------------------
// this function counts the occurences of each of a list of prevented sequences
// as they occur in S. 
void count_prevented(char *S, char **pv_list, int *pv_lengths, int num_pv, int *pv_occurrences)
{ 
  int i,j,slen, pvlen;
  slen = strlen(S);
  for (i = 0; i < num_pv; i++) {  // iterate over prevented sequences
    pv_occurrences[i] = 0;        // zero out from previous counting (if it happened)
    pvlen = strlen(pv_list[i]);     
    for (j=0; j<slen-(pvlen-1); j++) { // iterate over possible frames
      strncpy(prevented_s, S + j, pvlen);  // prevented_s is a global, maxPrevened size char array
      prevented_s[pvlen] = '\0'; // null terminate prevented_s
      if (compatible_sequence_cmp(pv_list[i],prevented_s)==0) // prevented string detected
	   {pv_occurrences[i]++;}         // add occurrence of ith sequence
    }
  }
}

// --------------------------------------------------------------------------------
int count_equivalence_classes(char *S, int *wc, int *eq)
{ 
  int i, n, nq=0;
  // count number of unique base equivalent classes
  n=strlen(S); for (i=0; i<n; i++)  if (eq[i]==i+1 && (wc[i]>i+1 || wc[i]==-1)) nq++; 
  return nq;
}
// --------------------------------------------------------------------------------

int set_auto_spurious_weights(char *S, int *wc, int *eq) // returns bmax value and sets score weights
{
  // n_i are strand lengths;  c are complexes, sets of strands
  // weight for length-k intra-strand bindings is beta^(k-k0) where k0 = log_4 (max_i n_i)  
  // weight for length-k inter-strand intra-complex bindings is beta^(k-k1) where k1 = log_4 (max_c sum_{i in c} n_i)
  // weight for length-k inter-complex bindings is beta^(k-k2) where k2 = log_4 (sum_i n_i)
  // the idea here being that k0, k1, k2 represent the total number of unique subsequences needed

  int i, n, nq=0, nk0=0, nk1=0, nk2=0, ni=0, nc=0;

  i=0; while (S[i]!=0) {
    if (S[i] != ' ') { nk2++; ni++; nc++; nk0=max(nk0,ni); nk1=max(nk1,nc); }
    if (S[i] == ' ') { ni=0; }
    if (S[i] == ' ' && i>0 && S[i-1] == ' ') { nc=0; }
    i++;
  }

  spurious_intraS = pow(spurious_beta,0.0-log2(1.0*nk0)/2);
  spurious_interS = pow(spurious_beta,0.0-log2(1.0*nk1)/2);
  spurious_interC = pow(spurious_beta,0.0-log2(1.0*nk2)/2);

  W_spurious = 1.0;
  W_verboten = 64.0 / nk2;  // logic: three-letter verboten sequences occur once/64, so "typical" score should be 1.0

  // count number of unique base equivalent classes
  n=strlen(S); for (i=0; i<n; i++)  if (eq[i]==i+1 && (wc[i]>i+1 || wc[i]==-1)) nq++; 
  if (!quiet) printf("Automatic: counted %d unique base equivalence classes.\n",nq);

  return (attempt_factor*3*nq); // for bmax so that optimizer 'gets bored' after 
                              // attempt_factor-fold coverage of all possible single mutations
                              // since mutate() only returns one of 3 non-identity mutations this is
                              // attempt_factor * 3 * nq.
}

// -------------------------------------------------------------------------------------
// scoring function from score_spurious.m in DNAdesign
double score_spurious(char *S, int *wc, int *eq)
{
  // prefactors = [5^(-4) 5^(-5) 5^(-6)];  % most weight to intramolecular
  // S0 = prefactors * (spurious(S, 3,8, wc,eq) * cumprod(5*ones(1,6))')  ;
  // S1 = prefactors * (spurious1(S, 5,10) * cumprod(5*ones(1,6))')  ;
  
  double factors[18];
  int *C; double score=0.0; int i;

  factors[0] =spurious_intraS; for (i=1; i<6; i++) factors[i]   = factors[i-1] *spurious_beta;
  factors[6] =spurious_interS; for (i=1; i<6; i++) factors[i+6] = factors[i+5] *spurious_beta;
  factors[12]=spurious_interC; for (i=1; i<6; i++) factors[i+12]= factors[i+11]*spurious_beta;

  C = spurious(0,S,3,8,wc,eq);
  for (i=0; i<18; i++) score+=factors[i]*C[i];

  C = spurious(1,S,5,10,NULL,NULL);
  for (i=0; i<18; i++) score+=factors[i]*C[i]*spurious_mm/spurious_beta/spurious_beta;

  return score;
}


// -------------------------------------------------------------------------------------
// scoring function from score_spurious.m in DNAdesign
// here S is over 1234, as is Swc
double score_spurious_reversible(char *S, char *Swc, int *wc, int *eq)
{
  
  int *C; double score=0.0; int i;

  // prefactors = [5^(-4) 5^(-5) 5^(-6)];  % most weight to intramolecular
  // S0 = prefactors * (spurious(S, 3,8, wc,eq) * cumprod(5*ones(1,6))')  ;
  // S1 = prefactors * (spurious1(S, 5,10) * cumprod(5*ones(1,6))')  ;

  C = spurious_reversible(0,S,Swc,3,8,wc,eq);
  for (i=0; i<18; i++) score+=factors_spurious[i]*C[i];

  C = spurious_reversible(1,S,Swc,5,10,NULL,NULL);
  for (i=0; i<18; i++) score+=factors_spurious[i]*C[i]*spurious_mm/spurious_beta/spurious_beta;

  return score;
}
// --------------------------------------------------------------------------

/*  wrapper to access  Vienna RNA package   c Ivo L Hofacker */

/*--------------------------------------------------------------------------*/

// pf:  0: mfe base pairs only
//      1: compute partition function base pair probabilities

char *fold_structure=NULL;  // keep fold data around until next use
char *fold_string=NULL;
int fold_length=0;

double foldit(char *S, int pf)
{
   int   i, length;  char *string, *structure;
   double energy=0.0, min_en=0.0;
   float kT, sfact=1.07;

   do_backtrack = 1;   // if partition function, get full matrix
                       // not just dG
   dangles=1;          // see man RNAfold... ignore restrictions 
                       // don't use dangles=1 with pf=1
   dangles=(pf?2:1);   // does this fix a bug?
   dangles=2;          // always....

   // temperature = 37; 

   length = strlen(S);

   if (length>fold_length) {
     if (fold_length>0) {
       free(fold_structure); free(fold_string);
       free_arrays(); free_pf_arrays();
     }
     // strange: I should be able to declare just 'length' here,
     // but if I do, then score_struct only works if 'quiet=TRUE' is *NOT* set
     // which, it turns out, has to do with executing the seq-NNNNN-seq
     // foldit call in score_nofold.   I don't understand, but using
     // '2*length' here seems to have fixed the problem.
     fold_string    = (char *) space(2*length+1); 
     fold_structure = (char *) space(2*length+1);
     initialize_fold(2*length);
     init_pf_fold(2*length);
     fold_length=2*length;
   }

   // printf("{%s}\n",S);

   string=fold_string; structure=fold_structure;
   strcpy(string,S);
   for (i = 0; i < length; i++) string[i] = toupper(string[i]);
   for (i = 0; i < length; i++) string[i] = (string[i]=='T')?'U':string[i];

   // pf needs min energy fold, for scale.
   min_en = fold(string, structure);  
   // printf("{[%s]}\n([%s]) min_en = %f\n",string,structure,min_en);
   // base_pairs[] now contains mfe structure base pairs
   if (pf) {
     kT = (temperature+273.15)*1.98717/1000.; /* in Kcal */
     pf_scale = exp(-(sfact*min_en)/kT/length);
     update_pf_params(length);
     energy = pf_fold(string, structure);
     if (debug) printf("\n[%s]\n[%s]\n energy = %f, mfe = %f\n",string,structure,energy,min_en); 

     // pr[] now contains base pairing probabilities 
   }

   return (pf?energy:min_en);
}

int num_strands(char *rS)                     // count number of strings
{
  char *S,*s,n=0; S=rS;
  while (*S != 0) {
   if ((s=strchr(S,' ')) != NULL) {
     if (s!=S) n=n+1;
     S = s+1;
   } else {
     n=n+1;
     S=rS+strlen(rS);
   }
  }
  return n;
}  

// -------------------------------------------------------------------------------------
double score_nofold(char *rS, int *wc, int *eq)
{
 char *S, *s, **si; int n=0; int i,j; char *string;
 double *dG; double mdG=0, mdGG=0, mddG=0, ndG=0, ndGG=0, nddG=0; double dg;

 string = (char *) malloc(sizeof(char)*2*strlen(rS)+10);

 n=num_strands(rS);

 si = (char **) malloc( sizeof(char *) * n);    // find string start indices
 S = rS; n=0;
 while (*S != 0) {
   if ((s=strchr(S,' ')) != NULL) {
     if (s!=S) { si[n]=S; n=n+1; }
     S = s+1;
   } else {
     si[n]=S; n=n+1;
     S=rS+strlen(rS);
   }
 }  
 for (i=1; i<n; i++) *(si[i]-1)=0;        // terminate strings

// for (i=0; i<n; i++) printf("%s\n",si[i]);


 dG = (double *) malloc( sizeof(double) * n);
 for (i=0; i<n; i++) if (nofold_include[i]) {
     dG[i] = foldit(si[i],0);      // mfe energy
     mdG = max(mdG, -dG[i]);  ndG += -dG[i];
 }  
 ndG=ndG/n;

// for (i=0; i<n; i++) printf("%6.4f ",-dG[i]); printf("\n\n");

 for (i=0; i<n; i++) if (nofold_include[i]) {
   for (j=0; j<n; j++) if (nofold_include[j]) {
    sprintf(string,"%sNNNNN%s",si[j],si[i]); 
    dg = foldit(string,0);      // mfe energy
    mdGG = max(mdGG, -dg); ndGG += -dg;
    mddG = max(mddG, -dg+dG[i]+dG[j]); nddG += -dg+dG[i]+dG[j];

//    printf("%6.4f ",-dg);
   }
//   printf("\n");
 }
 ndGG=ndGG/n/n; nddG=nddG/n/n;

 for (i=1; i<n; i++) *(si[i]-1)=' ';  // reconvert to single string

 free(si); free(string); free(dG);

// printf("[%f %f %f %f]\n",mdG, mdGf, mddG, V);

 return nofold_mdG*mdG + nofold_mdGG*mdGG + nofold_mddG*mddG +
        nofold_ndG*ndG + nofold_ndGG*ndGG + nofold_nddG*nddG;

}

// -------------------------------------------------------------------------------------
// with respect to intramolecular folds only,
// N - expected number of bases correctly paired  - # spaces in sequence 
//     == 0 if perfect
// no non-intramolecular interactions are considered at all.
double score_struct(char *rS, int *wc, int *eq)
{
 char *S, *s, **si; int n=0, N, Nk; int i,j,k,jjp; char *string;
 double expNbp=0.0, dG, score; 

 N = strlen(rS);

 string = (char *) malloc(sizeof(char)*2*N+10);

 S = rS;                                 // count number of strings
 while (*S != 0) {
   if ((s=strchr(S,' ')) != NULL) {
     if (s!=S) n=n+1;  S = s+1;
   } else {
     n=n+1;  S=rS+strlen(rS);
   }
 }  
 si = (char **) malloc( sizeof(char *) * n);    // find string start indices
 S = rS; n=0;
 while (*S != 0) {
   if ((s=strchr(S,' ')) != NULL) {
     if (s!=S) { si[n]=S; n=n+1; } S = s+1;
   } else {
     si[n]=S; n=n+1; S=rS+strlen(rS);
   }
 }  
 for (k=1; k<n; k++) *(si[k]-1)=0;        // terminate strings

 // k is strand #; i,j are within-strand indices 1..Nk; ii,jj are indices into S, 0..(N-1)
 // ii == (i-1)+(si[k]-rS);   i == (ii+1)-(si[k]-rS);  iip == ii+1;
 for (k=0; k<n; k++) { 
     Nk = strlen(si[k]);
     dG = foldit(si[k],1);      // pf energy; bp probs in pr[iindx[i]-j]]
     // pr is initially upper triangular, other side completely random
     for (i=1; i<=Nk; i++) for (j=1; j<i; j++) pr[iindx[i]-j]=pr[iindx[j]-i];
     for (i=1; i<=Nk; i++) {
       if ( (jjp=wc[(i-1)+(si[k]-rS)])==-1 ) {  // i should pair with nothing
         expNbp+=1.0;                          // counts as correct if it does
         for (j=1; j<=Nk; j++) expNbp-=pr[iindx[i]-j];
       } else {                                // i should pair with (eq class) jj       
         for (j=1; j<=Nk; j++) 
           if (eq[(j-1)+(si[k]-rS)] == jjp)     // j is in eq class jj
             expNbp+=pr[iindx[i]-j];            // counts as correct if paired
       }
       //       printf("%d->%d: %f\n",i,(j==Nk+1)?-1:j,expNbp);
     }
 }  

 for (k=1; k<n; k++) *(si[k]-1)=' ';  // reconvert to single string

 free(si); free(string);

 score = N - (n-1) - expNbp; // expected number of nt that *aren't* paired correctly.

 return score;

}

// -------------------------------------------------------------------------------------
// with respect to intramolecular folds only,
// prob( struct == target struct ) = exp(-E(struct,seq)/kT)/Z 
//     == 1.0 if perfect
// no non-intramolecular interactions are considered at all.
double score_prob(char *rS, int *wc, int *eq)
{
 char *S, *s, **si; int n=0, N, Nk; int i,j,k,jj,ii; char *string;
 double dG, dGtarget, score=0.0; double kT;

 N = strlen(rS);

 string = (char *) malloc(sizeof(char)*2*N+10);

 S = rS;                                 // count number of strings
 while (*S != 0) {
   if ((s=strchr(S,' ')) != NULL) {
     if (s!=S) n=n+1;  S = s+1;
   } else {
     n=n+1;  S=rS+strlen(rS);
   }
 }  
 si = (char **) malloc( sizeof(char *) * n);    // find string start indices
 S = rS; n=0;
 while (*S != 0) {
   if ((s=strchr(S,' ')) != NULL) {
     if (s!=S) { si[n]=S; n=n+1; } S = s+1;
   } else {
     si[n]=S; n=n+1; S=rS+strlen(rS);
   }
 }  
 for (k=1; k<n; k++) *(si[k]-1)=0;        // terminate strings

 // k is strand #; i,j are within-strand indices 1..Nk; ii,jj are indices into S, 0..(N-1)
 // ii == (i-1)+(si[k]-rS);   i == (ii+1)-(si[k]-rS); 
 for (k=0; k<n; k++) { 
     if (debug) printf("|%s|\n",si[k]);
     Nk = strlen(si[k]);
     dG = foldit(si[k],1);      // pf energy; bp probs in pr[iindx[i]-j]]
     // create *intramolecular* component of target structure
     //     if there are multiple wc partners, choose the lowest-indexed
     //     one within this strand.
     for (i=1; i<=Nk; i++) {
       ii=(i-1)+(si[k]-rS); jj=wc[ii]-1; j=(jj+1)-(si[k]-rS); 
       //       printf("i=%d ii=%d j=%d jj=%d\n",i,ii,j,jj);
       if      (jj!=-2 && 1<=j && j<i)    string[ii]=')';
       else if (jj!=-2 && i<j && j<=Nk)   string[ii]='(';
       else if (jj==-2)                   string[ii]='.';
       else if (j<1) { int j2, jj2;
           // i makes a bp to someone, see if anyone in strand is eq.
           // this can only be problematic for multi-strand designs.
                                          string[ii]='.';
          for (j2=1; j2<=Nk; j2++) {
	    jj2=(j2-1)+(si[k]-rS); 
            if      (eq[jj2]==jj && j2<i) string[ii]=')';
            else if (eq[jj2]==jj && i<j2) string[ii]='(';
	  }
       } else                             string[ii]='.';
     }
     string[Nk+(si[k]-rS)]=0;
     if (0 && debug) printf("|%s|\n",&string[si[k]-rS]);
     // ***** PROBLEM: wc does not always specify the desired secondary structure, 
     //       if eq has non-trivial values, or if pseuodknot is requested, etc.
     //       Thus we may have computed an impossible
     //       "target structure", and Vienna will crash.  EW & JMS 3/19/03
     // so: double check, & trim out bad things
     { int d; 
        for (d=0, i=1; i<=Nk; i++) {
          ii=(i-1)+(si[k]-rS);
          if (string[ii]=='(') d++;
          if (string[ii]==')') d--;
          if (d<0) { string[ii]='.'; d++; printf("#"); }
        }
        for (d=0, i=Nk; i>=1; i--) {
          ii=(i-1)+(si[k]-rS);
          if (string[ii]==')') d++;
          if (string[ii]=='(') d--;
          if (d<0) { string[ii]='.'; d++; printf("@"); }
        }
     }
     if (debug) printf("\n|%s|\n",&string[si[k]-rS]);

     dGtarget = energy_of_struct(si[k],&string[si[k]-rS]);
     kT = (temperature+273.15)*1.98717/1000.; /* in Kcal */
     if (debug) printf("[%s] T=%f, dGtarget=%f, dG=%f, p=%f\n",
            &string[si[k]-rS],temperature,dGtarget,dG,exp(-(dGtarget-dG)/kT));
     score += 1-exp(-(dGtarget-dG)/kT);
 }  

 for (k=1; k<n; k++) *(si[k]-1)=' ';  // reconvert to single string

 free(si); free(string);

 return score;

}

//------------------------------------------------------------------------------
// this is the score function invoked in the main optimization loop.
double score_all(char *rS, int *wc, int *eq, int num_tar, int **targeted_list)
{
  double V=0.0, Spu=0.0, Bnd=0.0, Nof=0.0, Str=0.0, Pro=0.0, Tar=0.0;

  if (W_verboten>0) V   = score_verboten(rS,wc,eq);
  if (W_spurious>0) Spu = score_spurious(rS,wc,eq);
  if (W_bonds>0)    Bnd = score_bonds(rS,wc,eq);
  if (W_nofold>0)   Nof = score_nofold(rS,wc,eq);
  if (W_struct>0)   Str = score_struct(rS,wc,eq);
  if (W_prob>0)     Pro = score_prob(rS,wc,eq);
  if (W_targeted>0)   Tar = score_targeted(rS,wc,eq, num_tar, targeted_list);

  return W_verboten*V + W_spurious*Spu + W_bonds*Bnd +
         W_nofold*Nof + W_struct*Str + W_prob*Pro + W_targeted*Tar;

}

//------------------------------------------------------------------------------
// this is the score function invoked in the main optimization loop.
// for use with spurious reversible.
// here rS is ACGT, but S and Sw are 1234
double score_all_reversible(char *rS, char *S, char *Sw, int *wc, int *eq, int num_tar, int **targeted_list)
{
  double V=0.0, Spu=0.0, Bnd=0.0, Nof=0.0, Str=0.0, Pro=0.0, Tar=0.0;

  if (W_verboten>0) V   = score_verboten(rS,wc,eq);
  if (W_spurious>0) Spu = score_spurious_reversible(S,Sw,wc,eq);
  if (W_bonds>0)    Bnd = score_bonds(rS,wc,eq);
  if (W_nofold>0)   Nof = score_nofold(rS,wc,eq);
  if (W_struct>0)   Str = score_struct(rS,wc,eq);
  if (W_prob>0)     Pro = score_prob(rS,wc,eq);
  if (W_targeted>0)   Tar = score_targeted_reversible(S,Sw,wc,eq, num_tar, targeted_list);

  return W_verboten*V + W_spurious*Spu + W_bonds*Bnd +
         W_nofold*Nof + W_struct*Str + W_prob*Pro + W_targeted*Tar;

}


// ----------------------------------------------------------------------------------------
// This function reports the score of all the relevant tests 
// as dictated by command line arguments---it is invoked once before 
// the optimization loop and once afterwards. It is a verbose accounting of 
// the components of the score function added up by 'score_all'

void test_evals(char *testS, int *testwc, int *testeq, int num_tar, int **targeted_list)
{
  int *C; int i,k; char *S,*s; double dG=0, V=0, Spu=0, Bnd=0, Nof=0, Str=0, Pro=0, Tar=0;
  int sh = (quiet==0 && watch==0);

  printf("\n");

  if (W_spurious>0) {

    C = spurious(0, testS, 3,8, testwc, testeq);
    if (spurious_equality) printf("spurious counts identity matches as well as WC matches.\n");
    else printf("spurious counts WC matches but not identity matches.\n");
    printf("spurious(testS, 3,8, testwc, testeq)\n");
    printf("C =  \n");
    for (i=0;i<3;i++) {
      for (k=3; k<=8; k++) {
	printf(" %5d", C[i*6 + k-3]);
      }
      printf("\n");
    }
    free(C);
    
    C = spurious(0, testS, 3,8, NULL, NULL);
    printf("spurious0(testS, 3,8)\n");
    printf("C =  \n");
    for (i=0;i<3;i++) {
      for (k=3; k<=8; k++) {
	printf(" %5d", C[i*6 + k-3]);
      }
      printf("\n");
    }
    free(C);
    
    
    C = spurious(1, testS, 5,10, NULL, NULL);
    printf("spurious1(testS, 5,10)\n");
    printf("C =  \n");
    for (i=0;i<3;i++) {
      for (k=5; k<=10; k++) {
	printf(" %5d", C[i*6 + k-5]);
      }
      printf("\n");
    }
    free(C);
  }

  printf("spurious: intraS = %6.4f, interS = %6.4f, interC = %6.4f\n", 
	 spurious_intraS, spurious_interS, spurious_interC);
  
  if (W_nofold>0 || W_struct>0 || W_prob>0) {
    int i, nofold_only = (W_struct==0 && W_prob==0);
    S = testS; i=0;
    while (*S != 0) {
      int last_one; last_one = ((s=strchr(S,' ')) == NULL);
      if (s!=S) {
        if (!last_one)  *s = 0; 
        // printf("%d: %s %d (%s)\n", i+1, S, nofold_include[i], last_one?"last":"more");
        if (nofold_only && nofold_include[i]) {
          dG = foldit(S,0);  printf("strand %d: %s: mfe = %6.4f; ",i+1,S,dG);
          dG = foldit(S,1);  printf("energy  = %6.4f\n",dG);
	}
        if (!last_one)  *s = ' '; 
        i++;
      }
      if (last_one) S=testS+strlen(testS); else S = s+1;
    }
  }

  // need to fill this in with a variety of targeted spurious tests.   
  if (W_targeted > 0) {
    printf("\n");
    printf("Targeted spurious score\n");

    C = targeted_spurious(0, testS, 3,8, testwc, testeq, num_tar,targeted_list);
    if (spurious_equality) printf("targeted spurious counts identity matches as well as WC matches.\n");
    else printf("targeted spurious counts WC matches but not identity matches.\n");
    printf("targeted spurious(testS, 3,8, testwc, testeq)\n");
    printf("C =  \n");
    for (i=0;i<2;i++) {
      for (k=3; k<=8; k++) {
	printf(" %5d", C[i*6 + k-3]);
      }
      printf("\n");
    }
    free(C);
    
    C = targeted_spurious(0, testS, 3,8, NULL, NULL,num_tar,targeted_list);
    printf("targeted spurious0(testS, 3,8)\n");
    printf("C =  \n");
    for (i=0;i<2;i++) {
      for (k=3; k<=8; k++) {
	printf(" %5d", C[i*6 + k-3]);
      }
      printf("\n");
    }
    free(C);
    
    C = targeted_spurious(1, testS, 5,10, NULL, NULL,num_tar,targeted_list);
    printf("targeted spurious1(testS, 5,10)\n");
    printf("C =  \n");
    for (i=0;i<2;i++) {
      for (k=5; k<=10; k++) {
	printf(" %5d", C[i*6 + k-5]);
      }
      printf("\n");
    }
    free(C);

  }


  if (sh || W_verboten>0) printf("** score_verboten score = %11.5f \n", V=score_verboten(testS,testwc,testeq) );
  if (sh || W_spurious>0) printf("** score_spurious score = %11.5f \n", Spu=score_spurious(testS,testwc,testeq) );
  if (sh || W_bonds>0) printf("** score_bonds    score = %11.5f \n", Bnd=score_bonds(testS,testwc,testeq) );
  if (sh || W_nofold>0) printf("** score_nofold   score = %11.5f \n", Nof=score_nofold(testS,testwc,testeq) );
  if (sh || W_struct>0) printf("** score_struct   score = %11.5f \n", Str=score_struct(testS,testwc,testeq) );
  if (sh || W_prob>0) printf("** score_prob     score = %11.5f \n", Pro=score_prob(testS,testwc,testeq) );
  if (sh || W_targeted>0) printf("** score_targeted score = %11.5f \n", Tar=score_targeted(testS,testwc,testeq,Ntar,targeted_list) );
  printf("** [verboten spurious bonds nofold struct prob targeted] = [%4.2f %4.2f %4.2f %4.2f, %4.2f, %4.2f, %4.2f]-weighted score = %11.5f \n", 
	 W_verboten, W_spurious, W_bonds, W_nofold, W_struct, W_prob, W_targeted,
	 W_verboten*V + W_spurious*Spu + W_bonds*Bnd + 
	 W_nofold*Nof + W_struct*Str + W_prob*Pro + W_targeted*Tar);
  
  printf("\n\n");

}


// ----------------------------------------------------------------------------------------
// Version of test_evals for reversible spurious
// This function reports the score of all the relevant tests 
// as dictated by command line arguments---it is invoked once before 
// the optimization loop and once afterwards. It is a verbose accounting of 
// the components of the score function added up by 'score_all'
// testS is over ACGT, but S1234, S1234wc are over 1234

void test_evals_reversible(char *testS, char *S1234, char *S1234wc, int *testwc, int *testeq, 
			   int num_tar, int **targeted_list)
{
  int *C; int i,k; char *S,*s; double dG=0, V=0, Spu=0, Bnd=0, Nof=0, Str=0, Pro=0, Tar=0;
  int sh = (quiet==0 && watch==0);

  printf("\n");

  if (W_spurious>0) {

    C = spurious_reversible(0, S1234, S1234wc, 3,8, testwc, testeq);
    if (spurious_equality) printf("spurious counts identity matches as well as WC matches.\n");
    else printf("spurious counts WC matches but not identity matches.\n");
    printf("spurious(testS, 3,8, testwc, testeq)\n");
    printf("C =  \n");
    for (i=0;i<3;i++) {
      for (k=3; k<=8; k++) {
	printf(" %5d", C[i*6 + k-3]);
      }
      printf("\n");
    }
    free(C);
    
    C = spurious_reversible(0, S1234, S1234wc,  3,8, NULL, NULL);
    printf("spurious0(testS, 3,8)\n");
    printf("C =  \n");
    for (i=0;i<3;i++) {
      for (k=3; k<=8; k++) {
	printf(" %5d", C[i*6 + k-3]);
      }
      printf("\n");
    }
    free(C);
    
    
    C = spurious_reversible(1, S1234, S1234wc, 5,10, NULL, NULL);
    printf("spurious1(testS, 5,10)\n");
    printf("C =  \n");
    for (i=0;i<3;i++) {
      for (k=5; k<=10; k++) {
	printf(" %5d", C[i*6 + k-5]);
      }
      printf("\n");
    }
    free(C);
  }

  printf("spurious: intraS = %6.4f, interS = %6.4f, interC = %6.4f\n", 
	 spurious_intraS, spurious_interS, spurious_interC);
  
  if (W_nofold>0 || W_struct>0 || W_prob>0) {
    int i, nofold_only = (W_struct==0 && W_prob==0);
    S = testS; i=0;
    while (*S != 0) {
      int last_one; last_one = ((s=strchr(S,' ')) == NULL);
      if (s!=S) {
        if (!last_one)  *s = 0; 
        // printf("%d: %s %d (%s)\n", i+1, S, nofold_include[i], last_one?"last":"more");
        if (nofold_only && nofold_include[i]) {
          dG = foldit(S,0);  printf("strand %d: %s: mfe = %6.4f; ",i+1,S,dG);
          dG = foldit(S,1);  printf("energy  = %6.4f\n",dG);
	}
        if (!last_one)  *s = ' '; 
        i++;
      }
      if (last_one) S=testS+strlen(testS); else S = s+1;
    }
  }
  
  if (W_targeted > 0) {
    printf("\n");
    printf("Targeted spurious score\n");

    C = targeted_spurious_reversible(0, S1234, S1234wc, 3,8, testwc, testeq, num_tar,targeted_list);
    if (spurious_equality) printf("targeted spurious counts identity matches as well as WC matches.\n");
    else printf("targeted spurious counts WC matches but not identity matches.\n");
    printf("targeted spurious(testS, 3,8, testwc, testeq)\n");
    printf("C =  \n");
    for (i=0;i<2;i++) {
      for (k=3; k<=8; k++) {
	printf(" %5d", C[i*6 + k-3]);
      }
      printf("\n");
    }
    free(C);
    
    C = targeted_spurious_reversible(0, S1234, S1234wc, 3,8, NULL, NULL,num_tar,targeted_list);
    printf("targeted spurious0(testS, 3,8)\n");
    printf("C =  \n");
    for (i=0;i<2;i++) {
      for (k=3; k<=8; k++) {
	printf(" %5d", C[i*6 + k-3]);
      }
      printf("\n");
    }
    free(C);
    
    C = targeted_spurious_reversible(1, S1234, S1234wc, 5,10, NULL, NULL,num_tar,targeted_list);
    printf("targeted spurious1(testS, 5,10)\n");
    printf("C =  \n");
    for (i=0;i<2;i++) {
      for (k=5; k<=10; k++) {
	printf(" %5d", C[i*6 + k-5]);
      }
      printf("\n");
    }
    free(C);


  }


  if (sh || W_verboten>0) printf("** score_verboten score = %11.5f \n", V=score_verboten(testS,testwc,testeq) );
  if (sh || W_spurious>0) printf("** score_spurious score = %11.5f \n", Spu=score_spurious_reversible(S1234,S1234wc,testwc,testeq) );
  if (sh || W_bonds>0) printf("** score_bonds    score = %11.5f \n", Bnd=score_bonds(testS,testwc,testeq) );
  if (sh || W_nofold>0) printf("** score_nofold   score = %11.5f \n", Nof=score_nofold(testS,testwc,testeq) );
  if (sh || W_struct>0) printf("** score_struct   score = %11.5f \n", Str=score_struct(testS,testwc,testeq) );
  if (sh || W_prob>0) printf("** score_prob     score = %11.5f \n", Pro=score_prob(testS,testwc,testeq) );
  if (sh || W_targeted>0) printf("** score_targeted score = %11.5f \n", Tar=score_targeted_reversible(S1234,S1234wc,testwc,testeq,Ntar,targeted_list) );
  printf("** [verboten spurious bonds nofold struct prob targeted] = [%4.2f %4.2f %4.2f %4.2f, %4.2f, %4.2f, %4.2f]-weighted score = %11.5f \n", 
	 W_verboten, W_spurious, W_bonds, W_nofold, W_struct, W_prob, W_targeted,
	 W_verboten*V + W_spurious*Spu + W_bonds*Bnd + 
	 W_nofold*Nof + W_struct*Str + W_prob*Pro + W_targeted*Tar);
  
  printf("\n\n");

}

// -------------------------------------------------------------------------

void load_input_files()
{ FILE *f; char c; int i, j, s1s, s1e, s2s, s2e,s1len,s2len,shift, stemp; double r; 
  int *targeted_pair_wc_old, *targeted_pair_eq_old;

//// Note: we strip off trailing ' ' from rS and St,
//// and delete corresponding entries in wc and eq
//// but YELL if any deleted entries are not -1 or 0 respectively.

 testS = malloc(maxN); // big enough?
 if (rS_filename!=NULL) {
  if ( (f = fopen(rS_filename,"r")) == NULL ) {
   fprintf(stderr,"Can't open rS file <%s>\n",rS_filename); exit(-1);
  } else {
   while ((c = fgetc(f)) != EOF) 
      if (index("ATCGatcg ",c)!=NULL) testS[NrS++]=c; 
   testS[NrS]=0;
   if (debug) printf("S = <%s>\n",testS);
   NrS = strlen(testS);
   while (testS[NrS-1]==' ' && NrS>0) { NrS--; testS[NrS]=0; }
   fclose(f);
  }
 }

 testSt = malloc(maxN);
 if (St_filename!=NULL) {
  if ( (f = fopen(St_filename,"r")) == NULL ) {
   fprintf(stderr,"Can't open St file <%s>\n",St_filename); exit(-1);
  } else {
   while ((c = fgetc(f)) != EOF) 
      if (index("ATCGatcgRYWSMKBDHVNrywsmkbdhvn ",c)!=NULL) testSt[NSt++]=c; 
   testSt[NSt]=0;
   if (debug) printf("St = <%s>\n",testSt);
   NSt = strlen(testSt);
   while (testSt[NSt-1]==' ' && NSt>0) { NSt--; testSt[NSt]=0; }
   fclose(f);
  }
 }

 testwc = calloc(maxN, sizeof(int));
 if (wc_filename!=NULL) {
  if ( (f = fopen(wc_filename,"r")) == NULL ) {
   fprintf(stderr,"Can't open wc file <%s>\n",wc_filename); exit(-1);
  } else {
   //   for (i=0; i<N; i++) { fscanf(f," %lf",&r); testwc[i]=r; }
   while (fscanf(f," %lf",&r)>0) testwc[Nwc++]=r; 
   if (debug) {
    printf("wc = <"); 
    for (i=0; i<Nwc; i++) { printf("%d ", testwc[i]); }
    printf("> \n"); 
   }
   for (i=0; i<Nwc; i++) if (testwc[i]==0 || testwc[i]<-1) 
     { fprintf(stderr,"error at position %d : value %d isn't allowed in wc\n", i+1, testwc[i]); exit(-1); }
   fclose(f);
  }
 }

 testeq = calloc(maxN, sizeof(int));
 if (eq_filename!=NULL) {
  if ( (f = fopen(eq_filename,"r")) == NULL ) {
   fprintf(stderr,"Can't open eq file <%s>\n",eq_filename); exit(-1);
  } else {
   //   for (i=0; i<N; i++) { fscanf(f,"%lf",&r); testeq[i]=r; }
   while (fscanf(f," %lf",&r)>0) testeq[Neq++]=r; 
   if (debug) {
    printf("eq = <"); 
    for (i=0; i<Neq; i++) { printf("%d ", testeq[i]); }
    printf("> \n"); 
   }
   while (testeq[Neq-1]==0 && Neq>0) Neq--;
   for (i=0; i<Neq; i++) if (testeq[i]<0) 
     { fprintf(stderr,"error at position %d : value %d isn't allowed in eq\n", i+1, testeq[i]); exit(-1); }
   fclose(f);
  }
 }

 // load in pv_filename, which should be a series of sequences separated by individual lines

 if (pv_filename!=NULL) {
   prevented_s = malloc(maxPrevented);
   if ( (f = fopen(pv_filename,"r")) == NULL ) {
     fprintf(stderr,"Can't open pv (prevent) file <%s>\n",pv_filename); exit(-1);
   } else {
     // first figure out how big it is
     while (fscanf(f,"%s\n", prevented_s)!=EOF) {
       Npv++;
     }
     fclose(f);
     f = fopen(pv_filename,"r");
     prevented_list = malloc(sizeof(char *) * Npv);
     prevented_lengths = malloc(sizeof(int *) * Npv);
     prevented_violations = malloc(sizeof(long int *) * Npv);
     prevented_occurrences_S = malloc(sizeof(int *) * Npv);
     prevented_occurrences_St = malloc(sizeof(int *) * Npv);
     j = 0;
     while (fscanf(f,"%s\n",prevented_s)!=EOF) {
       prevented_list[j] = malloc(sizeof(char) * (strlen(prevented_s) + 1));      
       strcpy(prevented_list[j],prevented_s);
       prevented_lengths[j] = strlen(prevented_list[j]);
       prevented_violations[j] = 0;
       prevented_occurrences_S[j] = 0;
       prevented_occurrences_St[j] = 0;
       make_upper(prevented_list[j]);  // convert to upper case if necessary
       printf("prevented string number %d: %s is %d bases long.\n", j+1, prevented_list[j], prevented_lengths[j]); 
       j++;
     }
     fclose(f);
   } // file opened
 } 


 // load in tar_filename 

 if (tar_filename!=NULL) {
   targeted_s1 = malloc(maxTargeted);
   targeted_s2 = malloc(maxTargeted);
   targeted_pair = malloc(maxTargeted);
   targeted_pair_comp = malloc(maxTargeted);
   targeted_s1t = malloc(maxTargeted);
   targeted_s2t = malloc(maxTargeted);
   targeted_pairt = malloc(maxTargeted);

   if ( (f = fopen(tar_filename,"r")) == NULL ) {
     fprintf(stderr,"Can't open tar file <%s> (extra important that these targeted regions are noncomplementary) \n", tar_filename); exit(-1);
   } else {
     // first figure out how big it is
     while (fscanf(f,"%i %i %i %i\n", &s1s,&s1e,&s2s,&s2e)!=EOF) {
       Ntar++;
     }
     fclose(f);
     f = fopen(tar_filename,"r");
     global_targeted_list = malloc(sizeof(int *) * Ntar);

     j = 0;
     while (fscanf(f,"%i %i %i %i\n", &s1s,&s1e,&s2s,&s2e)!=EOF) {
       // if the pairs are out of order, swap them. 
       // functions elsewhere assume that s1 occurs before s2
       if (s1s > s2s) {
	 stemp = s1s; s1s = s2s; s2s = stemp;
	 stemp = s1e; s1e = s2e; s2e = stemp;
       }
       global_targeted_list[j] = malloc(sizeof(int) * 4);      
       global_targeted_list[j][0] = s1s; global_targeted_list[j][1] = s1e; 
       global_targeted_list[j][2] = s2s; global_targeted_list[j][3] = s2e;
       printf("targeted spurious constraint %d is runs from %d to %d and %d to %d.\n",
	      j+1,global_targeted_list[j][0],global_targeted_list[j][1], 
	      global_targeted_list[j][2], global_targeted_list[j][3]); 
       if ((s1s>s1e) || (s2s>s1e)) {
         printf("Warning: target pair %d has a subsequence that ends before it starts.\n",j+1);
       }
       if ((s1s>=s2s) || (s1e>=s2s)) {
	 printf("Warning: overlapping subequences for targeted spurious pair %d.\n", j+1);}

       j++;
     }  // load in the target pairs
     fclose(f);
     printf("Warning: targeted spurious assumes all bases that pair do so uniquely---\n");
     printf("         it cannot handle wc files that specify a subsequence that is \n");
     printf("         complementary to two different subsequences.\n");



   } // we've got an open target file.
 } // a target file was specified


 // check length of wc compared to eq, sequence, and sequence template

 while (Nwc > MAX(Neq,MAX(NrS,NSt)) && MAX(Neq,MAX(NrS,NSt))>0) {
   if (testwc[Nwc-1]==-1) Nwc--;
   else {
     fprintf(stderr,"wc(%d) is longer than sequence(%d)/template(%d)/eq(%d), with non-empty entries!!\n",Nwc,NrS,NSt,Neq);
     exit(-1);
   }
 }

 N=MAX(N,MAX(MAX(NrS,NSt),MAX(Neq,Nwc)));
 if (N==0) {
   fprintf(stderr,"Zero-length sequence.  Aborting.  Try --help.  \n");
   exit(-1);
 }

 // set default array values for rS, St, wc, eq
 if (St_filename==NULL) {
   NSt=N; for (i=0; i<N; i++) testSt[i]='N'; testSt[N]=0;
 }
 if (rS_filename==NULL) {
   NrS=NSt; for (i=0; i<NrS; i++) testS[i]=randbasec(testSt[i]); testS[NrS]=0;
 }
 if (wc_filename==NULL) {
   Nwc=N; for (i=0; i<N; i++) testwc[i]=-1;
 }
 if (eq_filename==NULL) {
   Neq=N; for (i=0; i<N; i++) testeq[i]=(testS[i]==' ')?0:(i+1); 
 }


 if (N != NrS) 
   { fprintf(stderr,"rS is not max length, %d!!!\n",N); exit(-1); }
 if (N != NSt) 
   { fprintf(stderr,"St is not max length, %d!!!\n",N); exit(-1); }
 if (N != Nwc) 
   { fprintf(stderr,"wc is not max length, %d!!!\n",N); exit(-1); }
 if (N != Neq) 
   { fprintf(stderr,"eq is not max length, %d!!!\n",N); exit(-1); }

 // make corrections to defaults for ' ' separators
 for (i=0; i<N; i++) 
    if (testSt[i]==' ' || testS[i]==' ' || testeq[i]==0) 
      { testSt[i]=' '; testS[i]=' '; testeq[i]=0; testwc[i]=-1; }

 for (i=0; i<N; i++) if (testeq[i]==0 && testS[i]!=' ')
   { fprintf(stderr,"error at position %d : eq can be 0 only in space between strands\n",i+1); exit(-1); }

 // build up wc and eq tables for target pairs (these will not change over the course of the
 // optimization so they can be made now.       

 if (tar_filename!=NULL) {

   targeted_pair_wc = malloc(sizeof(int *) * Ntar);
   targeted_pair_eq = malloc(sizeof(int *) * Ntar);
   targeted_pair_wc_old =  malloc(maxTargeted*sizeof(int));
   targeted_pair_eq_old =  malloc(maxTargeted*sizeof(int));

   for (i = 0; i < Ntar; i++) {

     targeted_pair_wc[i] = malloc(maxTargeted*sizeof(int));
     targeted_pair_eq[i] = malloc(maxTargeted*sizeof(int));
     
     // get indices of pair i of targeted strings,
     // they start with 1 (from matlab) so  subtract 1 so they match C array convention.
     s1s = global_targeted_list[i][0]-1;
     s1e = global_targeted_list[i][1]-1;
     s2s = global_targeted_list[i][2]-1;
     s2e = global_targeted_list[i][3]-1;
     s1len = s1e - s1s + 1;
     s2len = s2e - s2s + 1;
     // copy targeted strings out of S, to targeted_s1, targeted_s2, global strings of size maxTargeted
     strncpy(targeted_s1, testS + s1s, s1len);
     strncpy(targeted_s2, testS + s2s, s2len);
     strncpy(targeted_s1t, testSt + s1s, s1len);
     strncpy(targeted_s2t, testSt + s2s, s2len);
     // null terminate
     targeted_s1[s1len] = '\0';
     targeted_s2[s2len] = '\0';
     targeted_s1t[s1len] = '\0';
     targeted_s2t[s2len] = '\0';
     // catenate to create a two strand system
     strcpy(targeted_pair, targeted_s1);
     strcat(targeted_pair, " ");
     strcat(targeted_pair, targeted_s2);
     strcpy(targeted_pairt, targeted_s1t);
     strcat(targeted_pairt, " ");
     strcat(targeted_pairt, targeted_s2t);
     if ((strlen(targeted_pair) != s1len + s2len + 1) || (strlen(targeted_pair) != strlen(targeted_pairt)))
       {printf("Warning: error catenating targets for targeted spurious.\n");}


     // To construct new wc and eq without mistakes, we do it in a two-stage
     // and somewhat unwieldy process, of first shifting both strings reference frame,
     // in unison, and second shifting the frame of S2 down so is next to s1. 
     // Also, this code assumes uniqueness of wc pairings
     // and so won't work in general for circuits, etc. 
     
     // So that we aren't terminally confused about the different between a matlab index and
     // a C index, we subtract 1 from the eq and wc at the same time, we'll put this back later.
     // We do this at the same time we execute the first step,
     // shifting wc and eq by the index of s1s, and writing into targeted_s1_wc, etc. 
     // All intervening sequence, etc, is just ignored.     
     // For debugging purpose we just shift the EQ and WC constraints for s1 and s2 down
     // completely unmodified into targeted_pair_eq_old and targeted_pair_wc_old. They
     // should have wc and eq constraints that line up in a sensible way with the modifed ones.

     for (j=0; j<s1len; j++) { // shift down s1 
       targeted_pair_wc[i][j] = testwc[s1s +j] - s1s - 1;  // j = s1s+j - s1s
       targeted_pair_wc_old[j] = testwc[s1s + j];
       targeted_pair_eq[i][j] = testeq[s1s +j] - s1s - 1;
       targeted_pair_eq_old[j] = testeq[s1s + j];
     }
     for (j=0; j<s2len; j++) { // shift down s2 by the same amount
       targeted_pair_wc[i][s2s + j - s1s] = testwc[s2s +j] - s1s - 1;
       targeted_pair_wc_old[s1len + 1 + j] = testwc[s2s +j];
       targeted_pair_eq[i][s2s + j - s1s] = testeq[s2s +j] - s1s - 1;
       targeted_pair_eq_old[s1len + 1 + j] = testeq[s2s +j];
     }
     // now:
     s2s = s2s - s1s;
     s2e = s2s - s1s;
     s1e = s1e - s1s;
     s1s = 0;
     
     // now s2 has the appropriate final wc, eq, pointers to s1, since s1 is in 
     // it's final resting place, and s1 has the appropriate final eq pointers since
     // they should all be <= positions in s1. So we only need to fix wc pointers from
     // s1 to s2 and shift them down as we shift s2 down to be next to s1.
     
     // shift down s2 to its final resting place.
     // or potentially move it up one position, to leave a place for a space.
     // s2 will run from s1e + 2, to s1e + 1 + slen2
     // currently at s2s. So pos - (s2s - s1e - 2). When pos = s2s it starts correctly.
     
     shift = s2s - s1e - 2;
     for (j=0; j<s1len; j++) { // first shift wc constraints in s1
       if ((targeted_pair_wc[i][j] >= s2s) && (targeted_pair_wc[i][j] <= s2e)) { // s1 points to s2
	 targeted_pair_wc[i][j] = targeted_pair_wc[i][j] - shift;
       }
       else if ((targeted_pair_wc[i][j] >= s1s) && (targeted_pair_wc[i][j] <= s1e)) { // s1 points to s1
	 // OK do nothing
       }
       else { targeted_pair_wc[i][j] = -2;}  // s1 points outside of s1 and s2
     }
     
     for (j=s2s; j<=s2e; j++) { // now safely shift down s2
       
       if ((targeted_pair_wc[i][j] >= s1s) && (targeted_pair_wc[i][j] <= s1s)) { // s2 points to s1
	 targeted_pair_wc[i][j - shift] = targeted_pair_wc[i][j]; }
       else if ((targeted_pair_wc[i][j] >= s2s) && (targeted_pair_wc[i][j] <= s2s)) { // s2 points to s2
	 targeted_pair_wc[i][j - shift] = targeted_pair_wc[i][j] - shift;
       }
       else {
	 targeted_pair_wc[i][j - shift] = -2; // s2 points outside of s1 and s2
       }
     }
     
     s2s = s2s - shift;
     s2e = s2e - shift;
     if ((s1len + 1 + s2len) != (s2e + 1)) 
       {printf("Warning: error in length of catenated target pair %d.\n",i);}
     
     // add space between s1 and s2.
     targeted_pair_wc[i][s1e+1] = -2; //will become -1
     targeted_pair_eq[i][s1e+1] = -1; // will become 0
     targeted_pair_wc_old[s1len] = -1;
     targeted_pair_eq_old[s1len] = 0;
     
     // now run through s1/s2 to see that they have eq and wc that fall in 
     // the appropriate ranges and get rid of any pointers to the outside.       
     // also shift all vaules up by 1 to get back to the wc, eq standard...
     for (j=0; j<=s2e; j++) {
       if (((targeted_pair_wc[i][j]>=0) && (targeted_pair_wc[i][j]<=s2e) && (targeted_pair_wc[i][j] != s1e)) 
	   || targeted_pair_wc[i][j] == -2 ) {
	 //that's OK, do nothing
       }
       else { 
	 printf("Warning, WC constraint for targeted pair %d is out of range after shifting.\n",i);
       }
       
       if (((targeted_pair_eq[i][j]>=0) && (targeted_pair_eq[i][j]<=j) && (targeted_pair_eq[i][j]!=s1e)) 
	   || ((targeted_pair_eq[i][j] == -1) && (j == s1e)) ) {
	 //that's OK, do nothing
       }
       else { 
	 printf("Warning, EQ constraint for targeted pair %d is out of range after shifting.\n",i);
       }    
       targeted_pair_wc[i][j]++;
       targeted_pair_eq[i][j]++;
     } // check EQ and WC are OK

     // print out the targeted pairs and their adjusted WC and EQs
     printf("temp,  target %d: ",i);     
     for (j=0; j<= s2e; j++) {
       printf("   %c",targeted_pairt[j]);
     }
     printf("\n");

     printf("ACGT,  target %d: ",i);     
     for (j=0; j<= s2e; j++) {
       printf("   %c",targeted_pair[j]);
     }
     printf("\n");

     printf("EQ,    target %d: ",i);     
     for (j=0; j<= s2e; j++) {
       printf("%d4",targeted_pair_eq[i][j]);
     }
     printf("\n");

     printf("EQold, target %d: ",i);     
     for (j=0; j<= s2e; j++) {
       printf("%d4",targeted_pair_eq_old[j]);
     }
     printf("\n");

     printf("WC,    target %d: ",i);     
     for (j=0; j<= s2e; j++) {
       printf("%d4",targeted_pair_wc[i][j]);
     }
     printf("\n");

     printf("WCold, target %d: ",i);     
     for (j=0; j<= s2e; j++) {
       printf("%d4",targeted_pair_wc_old[j]);
     }
     printf("\n");

   } // for i over all elements of the target list.

   free(targeted_pair_wc_old); free(targeted_pair_eq_old);

 } // if we loaded in a target list


 printf("%d p_rev size \n",N);
 p_rev = calloc(N,sizeof(int));  //% holds pointers to previous instance of subsequence at i
 S1234   = malloc(maxN);            // will hold testS in 1234 format.
 S1234wc = malloc(maxN);            // will hold WC of testS in 1234 format. 
 strcpy(S1234,testS);
 strcpy(S1234wc,testS);
 make1234(S1234,testS);
 WCstring1234(S1234,S1234wc);

}

// -------------------------------------------------------------------------
int main(int argc, char *argv[])
{
  int i; FILE *f;
  time_t t_start,t_now,tmax=0;  int smax=0; int bored, bmax=0; int automatic=0;
  int last_prevented_score=0, new_prevented_score;
  int mattempts = 0, current_mutations = 0, step_mutation_count=0;
  int fancy_args=1, trace=0; int DNAfold=1, RNAfold=0;
  int nofold_include_i = -1;  // if we ask for only specific strands included, which argv was that?

  
  P_rev = calloc(1<<(2*maxSubsequence),sizeof(int));  //% this is a large table which is used by spurious to 
                                                  // record which of 4^maxSubsequence instances of subsequence
                                                  // have been found in a sequence being analyzed.
  P_hist = calloc(maxN, sizeof(int)); // record of nonzero entries of P_rev
 // do we have fancy args or simple args?
 for (i=1;i<argc;i++) {
   if (index(argv[i],'=')==NULL) fancy_args=0;
 }

 if (!fancy_args && argc != 5) {
    fprintf(stderr,"Two command line formats are accepted for %s:\n",argv[0]);
    fprintf(stderr,"  %s rS_file St_file wc_file eq_file\n",argv[0]);
    fprintf(stderr,"or\n");
    fprintf(stderr,"  %s [option]*\n",argv[0]);
    fprintf(stderr,"where the options are:\n"
    " INPUT/OUTPUT\n"
    "  sequence=[file]      for initial (seed) nucleotide sequence (.rS file) [default: random]\n"
    "  template=[file]      for allowed nucleotides at each position (.St file) [default: all N]\n"
    "  wc=[file]            for Watson-Crick pairing constraints (.wc file) [default: all -1]\n"
    "  eq=[file]            for equality constraints (.eq file) [default: identity]\n"
    "  tar=[file]           for targeted spurious constraints (.tar file) [default: none]\n"
    "  pv=[file]            for prevented sequence strings (.pv file) [default: none]\n"
    "  output=[filename]    for final output sequence upon completion [default: stdout]\n"
    "  N=[value]            length of sequence, if no files are specified\n"
    "  quiet=TRUE           don't print anything, except errors & output\n"
    "  quiet=SCORES         output intial & final test scores, for all relevant score methods\n"
    "  quiet=WATCH          output running tally of best test scores, for chosen method\n"
    "  quiet=ALL            output all the above [default]\n"
    "  trace=ON             output the current sequence, at every improvement\n"
    " STOPPING CRITERIA\n"
    "  tmax=[value]         stop after specified time (in seconds) has elapsed\n"
    "  smax=[value]         stop after specified number of mutation iterations (steps) have been completed\n"
    "  bmax=[value]         stop after specified number of mutation iterations yield no change in score.\n"
    "  mmax=[value]         stop after specified number of attempts to find a mutation compatible with constraints.\n"
    "  attempt_factor=[value]  stop after specified multiplicity of possible mutations [default: 4X coverage]\n"
    " SCORING FUNCTION\n"
    "  score=automatic      set stopping criteria and spurious and verboten score weights based on design size\n"
    "  score=[verboten|spurious|bonds|nofold|struct|prob] what scoring function to use [default: spurious]\n"
    "     verboten: penalizes use of subsquences such as GGG, AAAA, WWWWWW, PuPyPuPyPuPy\n"
    "     spurious: penalizes pairs of undesired exact WC matches of lengths 3-8 (and 1-mm of lengths 5-10)\n"
    "     bonds:    prefers strong nearest-neighbor stacking energies for designed helices\n"
    "     nofold:   penalizes ssDNA folding dG, pairwise binding dG, and pairwise binding ddG\n"
    "     struct:   penalizes by the expected number of incorrectly-interacting nucleotides (ssDNA only)\n"
    "     prob:     prefers high ensemble probability of the target structure (ssDNA only)\n"
    "     setting 'score=xxx' is shorthand for 'W_xxx=1.0' and other weights 0.0; this can be further adjusted...\n"
    "  W_verboten=[value]               weight for verboten score [default: 0.0]\n"
    "  W_spurious=[value]               weight for spurious score [default: 1.0]\n"
    "  W_targeted=[value]               weight for targeted spurious score [default: 0.0]\n"
    "  W_bonds=[value]                  weight for bonds score [default: 0.0]\n"
    "  W_nofold=[value]                 weight for nofold score [default: 0.0]\n"
    "  W_struct=[value]                 weight for struct score [default: 0.0]\n"
    "  W_prob=[value]                   weight for prob score [default: 0.0]\n"
    "  spurious_beta=[value]            multiplicative penalty, per match base pair [default: 5.0]\n"
    "  spurious_intraS=[value]          initial penalty, intra-strand matches [default: 1/125]\n"
    "  spurious_interS=[value]          initial penalty, inter-strand matches [default: 1/625]\n"
    "  spurious_interC=[value]          initial penalty, inter-complex matches [default: 1/3125]\n"
    "  spurious_equality=[0|1]          penalize identity matches, not just WC complementary matches [default: 1]\n"
    "  reversible_spurious=[0|1]        uses a reversible version of spurious to save memory allocation [default: 1]\n"
    "  nofold_set=[id],[id],...,[id]    subset of strand id's for nofold calculations [default: all]\n"
    "  nofold_mdG=[value]               penalty for max single strand mfe [default: 1.0]\n"
    "  nofold_mdGG=[value]              penalty for max strand-pair mfe [default: 1.0]\n"
    "  nofold_mddG=[value]              penalty for max strand-pair relative mfe [default: 1.0]\n"
    "  nofold_ndG=[value]               penalty for mean single strand mfe [default: 1.0]\n"
    "  nofold_ndGG=[value]              penalty for mean strand-pair mfe [default: 1.0]\n"
    "  nofold_nddG=[value]              penalty for mean strand-pair relative mfe [default: 1.0]\n"
    "  param=[DNA|RNA]                  what Vienna RNAfold parameter set to use [default: DNA]\n"
    "  temperature=[value]              temperature for energy and partition-function calculations [default: 37]\n"
    "\n"
    "Notes:\n\n"
    "sequence is in { A C G T }*\n\n"
    "template is in { A C G T  R  Y  W  S  M  K  B   D   H   V   N   }*\n"
    "        meaning  A C G T AG CT AT CG AC GT CGT AGT ACT ACG ACGT \n\n"
    "        spaces indicate strand (' ') and complex ('  ') boundaries\n\n"
    "If the seed sequence does not match the template base types, it will be constrained to\n"
    "to match the template sequence and a warning will be generated.\n"
    "Lowercase input is converted to uppper case.\n\n"
    "eq and wc are lists of integers, indexes into the sequence string, starting at 1\n\n"
    "eq(i) is the min of all bases constrained\n     to have the same value as i (' ' bases have eq(i)==0)\n\n"
    "wc(i) is the min of all bases constrained\n     to be complementary to i    (or -1 if there are none)\n\n"
    "Strand id's are numbered starting at 1.\n\n"
    "a prevent list has one sequence per line, in {ACGTRYWSMKBDHVN}*\n"
    "\n"
    "tmax is used to set the stopping condition based on time in seconds.\n"
    "     Variables 't_start' and 't_now' track the time.\n"
    "smax is used to set the stopping condition based on the number of mutations tested.\n"
    "     Variable 'steps' is limited by smax.\n"      
    "bmax is used to set the stopping condition based on a change of score --\n"
    "     how many successful mutations can be made without a score change before stopping.\n"
    "     Variable 'bored' is limited by bmax.\n"
    "mmax is used to set the stopping condition based on whether mutations\n"
    "     compatible with the prevent list and wc/eq constraints can be found at all\n"
    "     by the function that generates mutations -- \n"
    "     how many unsuccessful attempts at finding a mutation can be made before stopping.\n"
    "     Variable 'mattempts' is limited by mmax.\n"
    "\n"
    "Regarding sequence optimization:\n"
    "If all four limits are defined, optimization is stopped if any _one_ of the four limits is met.\n"
    "Similarly, for any defined subset of these four limits, stopping occurs when any \n" 
    "member of the subset is first satisfied.\n"
    "Only mmax is defined by default. It must be defined so that the mutation generator \n"
    "does not enter an infinite loop. When the mutation generator can no longer \n"
    "generate any mutations for the optimizer to score, because of the constraints placed on it\n"
    "there is no reason for the optimizer to continue.\n"
    "mmax is set to (attempt_factor * Nequiv *3), where attempt_factor is the X-fold coverage\n"
    "of each of Neqiv equivalence classes of unique bases, which can have one of 3 different\n"
    "mutations from the current base. By default attempt_factor=4 for 4X coverage.\n"
    "\n" 
    "When the program is started without any limits in the argument list, \n"
    "the program runs, outputting the sequence S every 300 'steps' where each step\n"
    "involves the evaluation of the score on a valid mutation, returned by the mutation generator.\n"
    "Optimization continues until (potentially) the mutation generator fails to make single base\n" 
    "changes that satisfy the constraints for mmax tries. Optimization will continue indefinitely\n"
    "if it is easy to make valid mutations (even if they make the score worse; the change to the\n"
    "sequence will just be rejected). The step at which the last improvement was made is reported\n"
    "so that the user can more easily identify putative local minima.\n"
    "\n"
    "Unless a prevent list complicates the constraints on a seqence, mmax is not usually the\n"
    "limit that causes the optimization loop to terminate. Which factor _is_ causing termination\n"
    "can be determined by comparing the relevant variable to its limit for each step at which\n"
    "a valid mutation (neutral or improving) is reported. If bored approaches or equals bmax\n"
    "it is likely to be responsible for termination given the input parameters (in similar runs).\n"
    "If mattempts approaches or equals mmax, then it is likely to be responsible; under this \n"
    "condition, little optimization occurs, and it may be useful to relax prevent constraints \n"
    "or sequence constraints in the template.\n"
    "\n"
    "Regarding prevent lists:\n"
    "When a prevent list (.pv file) is defined, a separate pre-processing phase may be invoked \n"
    "to remove prevented sequences from the seed sequence. As for the optimizer, this is done \n"
    "by random mutation. The template sequence (.St) may force a set of prevented sequences to \n"
    "occur in particular locations -- this is so that a restriction site, for example, \n"
    "may be specified to occur in one and only one place (by inclusion in the template sequence\n"
    "and on the prevent list). In general, the template sequence may force the occurrence of additional\n"
    "prevented sequences based on its structure. Thus the pre-processing loop seeks to minimize \n"
    "the 'prevent score' (the difference between the minimum number of prevented sequences _forced_\n"
    "by the template sequence and the number that occur in the seed.)\n"
    "During this pre-processing stage, mmax is also used as a limit on the number of invalid\n"
    "mutations that can be made before the pre-processing stage stops. Note that the mutation\n"
    "generator respects the prevent list, sequence template, wc and eq lists and will not \n"
    "generate a mutation that violates them.\n"
    "\n"
    "Regarding scoring functions:\n"
    "score = automatic sets bmax similarly to mmax as (attempt_factor * Neqiv *3) which is\n"
    "    (12 * Nequiv) unless 'attempt_factor' is specified.\n"
    "See DNAdesign score_verboten.m, score_spurious.m and score_nofold.m for scoring function definitions\n"
    "     'verboten', 'spurious' and 'nofold'.\n"
    "'struct' and 'prob' use the RNAfold partition function for single strands only.\n"
    "'nofold' uses the RNAfold mfe, and approximates bimolecular binding energies by adding a hairpin linker.\n"
    "\n"
    "Weighted combinations of scoring functions can be specified using the 'W_*' parameters.\n"
    "\n"
    "General:\n"
    "   It is useful to run the program multiple times with the same parameters to get a feel for\n"
    "   how the constraints effect the extent of the optimization, and which parameters are\n"
    "   are the dominant contributors to the score, given the constraints. Also, it is helpful\n"
    "   to correlate the 'quality' of sequence (however one defines it) across multiple runs\n"
    "   with the stopping parameters and constraints. Average base composition, especially, \n"
    "   (and obviously) has great effect on spurious and verboten scores. Change it (crudely)\n"
    "   by restricting the alphabet in the sequence template.\n"
    "\n"
    "   ********  To run the program just once with default parameters  **********\n"
    "   is probably to miss a much better sequence (for whatever specific purpose).\n"
    "\n"
    "Backwards incompatibilty:\n"
    "   smax used to be called imax, and so command lines invoking imax will break.\n");
    exit(-1);
 }

 init_rand();

 make_bad6mer();

 if (fancy_args) {
   // try to parse fancy args

   for (i=1;i<argc;i++) {
     if      (strncmp(argv[i],"sequence=",9)==0) rS_filename=&argv[i][9]; 
     else if (strncmp(argv[i],"template=",9)==0) St_filename=&argv[i][9]; 
     else if (strncmp(argv[i],"wc=",3)==0) wc_filename=&argv[i][3]; 
     else if (strncmp(argv[i],"eq=",3)==0) eq_filename=&argv[i][3]; 
     else if (strncmp(argv[i],"tar=",4)==0) tar_filename=&argv[i][4];
     else if (strncmp(argv[i],"pv=",3)==0) pv_filename=&argv[i][3];
     else if (strncmp(argv[i],"output=",7)==0) S_filename=&argv[i][7]; 
     else if (strncmp(argv[i],"N=",2)==0) N=atoi(&argv[i][2]);
     else if (strncmp(argv[i],"tmax=",5)==0) tmax=atoi(&argv[i][5]);
     else if (strncmp(argv[i],"smax=",5)==0) smax=atoi(&argv[i][5]);
     else if (strncmp(argv[i],"bmax=",5)==0) bmax=atoi(&argv[i][5]);
     else if (strncmp(argv[i],"mmax=",5)==0) mmax=atoi(&argv[i][5]);
     else if (strncmp(argv[i],"attempt_factor=",15)==0) attempt_factor=atoi(&argv[i][15]);
     else if (strncmp(argv[i],"score=automatic",15)==0) automatic=1;
     else if (strncmp(argv[i],"score=verboten",14)==0) { W_spurious=0.0; W_verboten=1.0; W_nofold=0.0; W_struct=0.0; W_prob=0.0; W_bonds=0.0; }
     else if (strncmp(argv[i],"score=spurious",14)==0) { W_spurious=1.0; W_verboten=0.0; W_nofold=0.0; W_struct=0.0; W_prob=0.0; W_bonds=0.0; }
     else if (strncmp(argv[i],"score=bonds",11)==0) { W_spurious=0.0; W_verboten=0.0; W_nofold=0.0; W_struct=0.0; W_prob=0.0; W_bonds=1.0; }
     else if (strncmp(argv[i],"score=nofold",12)==0) { W_spurious=0.0; W_verboten=0.0; W_nofold=1.0; W_struct=0.0; W_prob=0.0; W_bonds=0.0; }
     else if (strncmp(argv[i],"score=struct",12)==0) { W_spurious=0.0; W_verboten=0.0; W_nofold=0.0; W_struct=1.0; W_prob=0.0; W_bonds=0.0; }
     else if (strncmp(argv[i],"score=prob",10)==0) { W_spurious=0.0; W_verboten=0.0; W_nofold=0.0; W_struct=0.0; W_prob=1.0; W_bonds=0.0; }
     else if (strncmp(argv[i],"param=DNA",9)==0) { DNAfold=1; RNAfold=0; }
     else if (strncmp(argv[i],"param=RNA",9)==0) { RNAfold=1; DNAfold=0; }
     else if (strncmp(argv[i],"quiet=TRUE",10)==0) { quiet=1; watch=0; }
     else if (strncmp(argv[i],"quiet=SCORES",12)==0) { quiet=0; watch=0; }
     else if (strncmp(argv[i],"quiet=WATCH",11)==0) { quiet=1; watch=1; }
     else if (strncmp(argv[i],"quiet=ALL",9)==0) { quiet=0; watch=1; }
     else if (strncmp(argv[i],"trace=ON",8)==0) { trace=1;  }
     else if (strncmp(argv[i],"debug=TRUE",10)==0) debug=1;
     else if (strncmp(argv[i],"W_verboten=",11)==0) W_verboten=atof(&argv[i][11]);
     else if (strncmp(argv[i],"W_spurious=",11)==0) W_spurious=atof(&argv[i][11]);
     else if (strncmp(argv[i],"W_targeted=",11)==0) W_targeted=atof(&argv[i][11]);
     else if (strncmp(argv[i],"W_bonds=",8)==0) W_bonds=atof(&argv[i][8]);
     else if (strncmp(argv[i],"W_nofold=",9)==0) W_nofold=atof(&argv[i][9]);
     else if (strncmp(argv[i],"W_struct=",9)==0) W_struct=atof(&argv[i][9]);
     else if (strncmp(argv[i],"W_prob=",7)==0) W_prob=atof(&argv[i][7]);
     else if (strncmp(argv[i],"spurious_beta=",14)==0) spurious_beta=atof(&argv[i][14]);
     else if (strncmp(argv[i],"spurious_intraS=",16)==0) spurious_intraS=atof(&argv[i][16]);
     else if (strncmp(argv[i],"spurious_interS=",16)==0) spurious_interS=atof(&argv[i][16]);
     else if (strncmp(argv[i],"spurious_interC=",16)==0) spurious_interC=atof(&argv[i][16]);
     else if (strncmp(argv[i],"targeted_beta=",14)==0) targeted_beta=atof(&argv[i][14]);
     else if (strncmp(argv[i],"targeted_intraS=",16)==0) targeted_intraS=atof(&argv[i][16]);
     else if (strncmp(argv[i],"targeted_interS=",16)==0) targeted_interS=atof(&argv[i][16]);
     else if (strncmp(argv[i],"spurious_equality=",18)==0) spurious_equality=(atoi(&argv[i][18])>0);
     else if (strncmp(argv[i],"reversible_spurious=",20)==0) reversible_spurious=(atoi(&argv[i][20])>0);
     else if (strncmp(argv[i],"temperature=",12)==0) temperature=atof(&argv[i][12]);
     else if (strncmp(argv[i],"nofold_include=",15)==0) nofold_include_i=i;
     else if (strncmp(argv[i],"nofold_mdG=",11)==0) nofold_mdG=atof(&argv[i][11]);
     else if (strncmp(argv[i],"nofold_mdGG=",12)==0) nofold_mdGG=atof(&argv[i][12]);
     else if (strncmp(argv[i],"nofold_mddG=",12)==0) nofold_mddG=atof(&argv[i][12]);
     else if (strncmp(argv[i],"nofold_ndG=",11)==0) nofold_ndG=atof(&argv[i][11]);
     else if (strncmp(argv[i],"nofold_ndGG=",12)==0) nofold_ndGG=atof(&argv[i][12]);
     else if (strncmp(argv[i],"nofold_nddG=",12)==0) nofold_nddG=atof(&argv[i][12]);
     else { 
       fprintf(stderr,"Cannot parse argument <%s>.  Aborting. \n\n",argv[i]); exit(-1); 
     }
   }
 
 } else { 
   // otherwise, try to read all four input files
   rS_filename=argv[1]; St_filename=argv[2]; 
   wc_filename=argv[3]; eq_filename=argv[4];
 }

 // initialize RNAfold
 if ((f=fopen(VIENNA_DIR "default.par","r"))==NULL) { 
   fprintf(stderr, "Cannot open RNA param file at %s.  Recompile?\n", 
           VIENNA_DIR "default.par");
   exit(-1);
 } else { fclose(f); }
 if ((f=fopen(VIENNA_DIR "dna.par","r"))==NULL) { 
   fprintf(stderr, "Cannot open DNA param file at %s.  Recompile?\n", 
           VIENNA_DIR "dna.par");
   exit(-1);
 } else { fclose(f); }
 if (RNAfold) read_parameter_file(VIENNA_DIR "default.par");
 if (DNAfold) read_parameter_file(VIENNA_DIR "dna.par");

 // load in the files!
 load_input_files();

 //////////////////////////////////////////////////////////////
 // information was read in.  now do the designing!

 // condition the input to uppercase
 make_upper(testS); make_upper(testSt);

 // set up weights and factors....
 if (automatic) bmax=set_auto_spurious_weights(testS,testwc,testeq);

 //compute weights for spurious.
 // prefactors = [5^(-4) 5^(-5) 5^(-6)];  % most weight to intramolecular
 // S0 = prefactors * (spurious(S, 3,8, wc,eq) * cumprod(5*ones(1,6))')  ;
 // S1 = prefactors * (spurious1(S, 5,10) * cumprod(5*ones(1,6))')  ;
 factors_spurious[0] =spurious_intraS; for (i=1; i<6; i++) factors_spurious[i]   = factors_spurious[i-1] *spurious_beta;
 factors_spurious[6] =spurious_interS; for (i=1; i<6; i++) factors_spurious[i+6] = factors_spurious[i+5] *spurious_beta;
 factors_spurious[12]=spurious_interC; for (i=1; i<6; i++) factors_spurious[i+12]= factors_spurious[i+11]*spurious_beta;

 factors_targeted[0] =targeted_intraS; for (i=1; i<6; i++) factors_targeted[i]   = factors_targeted[i-1] *targeted_beta;
 factors_targeted[6] =targeted_interS; for (i=1; i<6; i++) factors_targeted[i+6] = factors_targeted[i+5] *targeted_beta;

 { int ns, i, j, s;  char *p;
   ns = num_strands(testS);
   nofold_include = (char *) malloc( ns * sizeof(char) );
   // by default, if nofold scoring is used, apply it to all strands
   // but if specific strands are asked for, default the others to not be used
   for (j=0; j<ns; j++) nofold_include[j]=(nofold_include_i == -1);
   // printf("nofold: parsing %s\n", argv[nofold_include_i]);
   if (nofold_include_i != -1) {  // find the strands specified to include
     i = nofold_include_i;  j=0;
     p = &argv[i][15];
     while (index(p,',') != NULL) {
       s = atoi(p); nofold_include[s-1] = 1;  p = index(p,',')+1;
       // printf("nofold: include strand %d\n", s);
     }
     s = atoi(p); nofold_include[s-1] = 1; // get last one
   }
 }

 printf("Attempt factor for mutations is %dX coverage of each equivalence class.", attempt_factor);

 i = template_match(testS,testSt);     
 if (i != strlen(testS))
   printf("\nWarning: seed sequence didn't match template at %d places and had to be matched to the template.\n", (strlen(testS)-i));

 Nequiv = count_equivalence_classes(testSt, testwc, testeq);
 Nonspace = count_nonspace(testSt);
 mmax = attempt_factor*3*Nequiv;

 constrain(testS,testwc,testeq);
 if (!quiet) printf("\nconstrained S = <\n%s\n>  N=%d \n",testS,N);
 // assumption is that the wc and eq constraints do not conflict with the sequence template!

 //--------------- Prevented sequence remover--------------------------------------------------------
 if (Npv != 0) { // remove all prevented sequences possible, count to see what is left due to the template
   count_prevented(testS, prevented_list, prevented_lengths, Npv, prevented_occurrences_S);   
   count_prevented(testSt, prevented_list, prevented_lengths, Npv, prevented_occurrences_St);   
   last_prevented_score = prevent_cmp(prevented_occurrences_S, prevented_occurrences_St, Npv);
   step_mutation_count = 0;
   if (!quiet)  {
     printf("Prevent score is the difference between the number of prevented sequences in\n");
     printf("the template St and the seed sequence S. (pre-optimized)\n");
     printf("Prevent score: %d (Initially)\n",last_prevented_score);}
   // do random mutations until the seed sequence has the same prevented sequences as the sequence template
   // or tired of mutations.
   while ((last_prevented_score !=0) && (mattempts < mmax)) {
     current_mutations = mutate(testS,testSt,testwc,testeq,prevented_list,prevented_lengths,Npv,prevented_violations);
     mattempts += current_mutations;
     step_mutation_count += current_mutations;
     count_prevented(testS, prevented_list, prevented_lengths, Npv, prevented_occurrences_S);   
     count_prevented(testSt, prevented_list, prevented_lengths, Npv, prevented_occurrences_St);   
     new_prevented_score = prevent_cmp(prevented_occurrences_S, prevented_occurrences_St, Npv);
     if ((!quiet) && (new_prevented_score < last_prevented_score)) {
          printf("Prevent score: %d   mutation attempts this step: %d   cumulative mutations: %d   mmax=%d\n",
		 new_prevented_score,step_mutation_count,mattempts,mmax);
	  step_mutation_count = 0;}
     last_prevented_score = new_prevented_score;
   }
   if (!quiet) printf("Prevent score: %d Final before optimization\n", last_prevented_score);
   if (!quiet) printf("\n prevent constrained S = <\n%s\n> \n",testS);
   for (i=0; i< Npv; i++) {
     if (prevented_occurrences_St[i] != 0) 
       if (!quiet) printf("Warning: Prevented sequence '%s' occurs %d time(s) due to the template. OK?\n", prevented_list[i], prevented_occurrences_St[i]);}
   if (prevent_cmp(prevented_occurrences_S, prevented_occurrences_St, Npv) !=0) {
     for (i=0; i<Npv; i++) { // what extra prevented sequences remain?
       if (prevented_occurrences_S[i] != prevented_occurrences_St[i]) 
	 printf("Warning: Number of prevented s '%s' in S (%d) doesn't match minimum possible, as dictated by template St (%d)\n", 
		prevented_list[i], prevented_occurrences_S[i], prevented_occurrences_St[i]);}
     printf("Warning: Gave up mutating seed sequence to purge prevented sequences after _total_ > %d mutations.\n",mmax);
     printf("         Look for warnings detailing extra prevented sequences after optimization.\n");   
     printf("         It may be by chance, or it may be logically impossible to purge more prevented sequences.\n");   
     printf("         If sequence examination determines it is by chance, increase 'attempt_factor' from default of 4.\n");
   }
   // clear out prevented_violations so that they can be counted during the main optimization loop.
   for (i=0; i<Npv; i++) {prevented_violations[i] = 0;}
 }  // end prevented sequence removal


 // -------------------------Main optimization code-----------------------------------------------
 if (!quiet) {
   if (reversible_spurious) {
     printf("Using a thermodynamically reversible version of spurious ;)\n");
     make1234(S1234,testS);
     WCstring1234(S1234,S1234wc);
     test_evals_reversible(testS,S1234,S1234wc,testwc,testeq,Ntar,global_targeted_list);
   }
   else {
     test_evals(testS,testwc,testeq,Ntar,global_targeted_list);
   }
 }
 if (count_nonspace(testSt) == count_fixed(testSt))
   { printf("All bases in sequence template fixed, no optimization performed.\n");}
 else {            // optimization block
   char *oldS; double old_score, new_score; int steps=0, last_improvement=0; bored=0;
  oldS = malloc(N); 
  mattempts = 0;
  if (reversible_spurious) {
    old_score = score_all_reversible(testS,S1234,S1234wc,testwc,testeq,Ntar,global_targeted_list); 
  }
  else {
    old_score = score_all(testS,testwc,testeq,Ntar,global_targeted_list); 
  }
  if (watch) printf("%8d steps, %8d seconds : score = %18.10f",0,0,old_score);
  if (watch && (smax==0) && (tmax==0) && (bmax==0))
    printf(" (mattempts=%d,mmax=%d)", mattempts,mmax);
  if (watch && (smax>0)) printf(" (smax=%d,mattempts=%d,mmax=%d)", smax,mattempts, mmax);
  if (watch && (tmax>0)) printf(" (tmax=%d,mattempts=%d,mmax=%d)", (int) tmax,mattempts, mmax);
  if (watch && (bmax>0)) printf(" (bored=%d,bmax=%d,mattempts=%d,mmax=%d)", bored, bmax,mattempts, mmax);
  if (watch) printf("\n");
  time(&t_start);  time(&t_now);
  while ( (tmax==0 || t_start+tmax>=t_now) && (smax==0 || steps<smax) && (bmax==0 || bored<bmax) && (mattempts<mmax)) { 
   steps++; 
   for (i=0; i<N; i++) oldS[i]=testS[i];

   if (reversible_spurious) {
     mattempts = mutate(testS,testSt,testwc,testeq,prevented_list,prevented_lengths,Npv,prevented_violations);
     make1234(S1234,testS);
     WCstring1234(S1234,S1234wc);
     new_score = score_all_reversible(testS,S1234,S1234wc,testwc,testeq,Ntar,global_targeted_list); 
   }
   else {
     mattempts = mutate(testS,testSt,testwc,testeq,prevented_list,prevented_lengths,Npv,prevented_violations);
     new_score = score_all(testS,testwc,testeq,Ntar,global_targeted_list); 
   }

   if (new_score <= old_score) {
     if (watch) printf("%8d steps, %8d seconds : score = %18.10f", 
		       steps, (int)(t_now-t_start), new_score);
     if (watch && (smax==0) && (tmax==0) && (bmax==0))
       printf(" (mattempts=%d,mmax=%d)", mattempts,mmax);
     if (watch && (smax>0)) printf(" (smax=%d,mattempts=%d,mmax=%d)", smax, mattempts, mmax);
     if (watch && (tmax>0)) printf(" (tmax=%d,mattempts=%d,mmax=%d)", (int) tmax, mattempts, mmax);
     if (watch && (bmax>0)) printf(" (bored=%d,bmax=%d,mattempts=%d,mmax=%d)", bored, bmax, mattempts, mmax);
     if (watch) printf("\n");
     if (trace) printf("%s\n",testS);
     if (new_score < old_score) {
       bored=0;
       last_improvement = steps;
       if (!(quiet && !watch) && !trace && steps%300 == 0 && tmax==0 && smax==0 && bmax==0) {
	 printf("\n%8d steps, %d seconds : score = %18.10f last improvement: this step.\n%s\n\n",steps,(int)(t_now-t_start), old_score,testS);
       }
     } 
     else {
       if (!(quiet && !watch) && !trace && steps%300 == 0 && tmax==0 && smax==0 && bmax==0) {
	 printf("\n%8d steps, %d seconds : score = %18.10f last improvement at step: %8d\n%s\n\n",steps,(int)(t_now-t_start), old_score,last_improvement,testS);
       }
     }
     old_score = new_score;
   } else {
     for (i=0; i<N; i++) testS[i]=oldS[i]; bored++; 
     if (!(quiet && !watch) && !trace && steps%300 == 0 && tmax==0 && smax==0 && bmax==0) {
       printf("\n%8d steps, %d seconds : score = %18.10f  last improvement at step: %8d\n%s\n\n",
	      steps,(int)(t_now-t_start), old_score, last_improvement, testS);
     }    
   }
   time(&t_now);
  }  
  if (watch) printf("%8d steps, %8d seconds : score = %18.10f FINAL\n", 
               steps-1, (int)(t_now-t_start), old_score);
  if (watch && (Npv !=0)) {
    printf("\n");
    for (i = 0; i < Npv; i++) printf("Prevented sequence '%s' had %ld violations during mutation.\n", prevented_list[i], prevented_violations[i]);
  }
 }  // end optimization block

 // -----------------------Prevented check post-optimization-------------------------------------------
 // check to see if the occurrences of prevented sequences in testS matches template after optimization
 if (!quiet && (Npv != 0)) { 
   count_prevented(testS, prevented_list, prevented_lengths, Npv, prevented_occurrences_S);   
   new_prevented_score = prevent_cmp(prevented_occurrences_S, prevented_occurrences_St, Npv);
   if (new_prevented_score < last_prevented_score) 
     printf("Prevent score improved during optimization from %d to %d.\n", last_prevented_score, new_prevented_score);
   else if (new_prevented_score == last_prevented_score) 
     printf("Prevent score remained the same during optimization: %d.\n", last_prevented_score);
   else printf("Warning: bad error, somehow the optimizer added a prevented sequence.\n");
   for (i=0; i<Npv; i++) {
     if (prevented_occurrences_S[i] != prevented_occurrences_St[i]) 
	 printf("Warning: Number of prevented s '%s' in S (%d) doesn't match minimum possible, as dictated by template St (%d)\n", 
		prevented_list[i], prevented_occurrences_S[i], prevented_occurrences_St[i]);}
 }

 // ***** print final sequence to output file *****

 if (!quiet) {
   if (reversible_spurious) {
     make1234(S1234,testS);
     WCstring1234(S1234,S1234wc);
     test_evals_reversible(testS,S1234,S1234wc,testwc,testeq,Ntar,global_targeted_list);
   }
   else {
     test_evals(testS,testwc,testeq,Ntar,global_targeted_list);
   }
 }

 if (S_filename==NULL)
   printf("%s\n",testS);
 else {
   if ( (f = fopen(S_filename,"w")) == NULL ) {
     fprintf(stderr,"Can't open output file <%s>\n",S_filename); exit(-1);
   } else {
     fprintf(f,"%s\n",testS);
     fclose(f);
   }
 }

 free(nofold_include);

 return 0;

}
