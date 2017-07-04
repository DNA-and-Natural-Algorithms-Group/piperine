/*

  Converted to C from Matlab by hand, EW, 4/2000 under name spuriousC.c

  command-line args enhanced, EW, 8/2002    

  Merged with spuriousCfold.c, EW, 8/2002.  Now *requires* the ViennaRNA package! 
  
  Removed the ViennaRNA dependencies, renamed as spuriousSSM.c, EW & CGE, 2/2015
  Also made mutations more efficient and reduced malloc/free load, resulting in 10x speed-up.

  gcc -Wall -O3 spuriousSSM.c -o spuriousSSM -lm

  Simple examples:

  ./spuriousSSM score=automatic template=examples/strands.st wc=examples/strands.wc eq=examples/strands.eq quiet=ALL
  ./spuriousSSM score=automatic template=examples/trans.st wc=examples/trans.wc eq=examples/trans.eq quiet=ALL
  ./spuriousSSM score=automatic template=examples/DAO.st wc=examples/DAO.wc eq=examples/DAO.eq quiet=ALL
  ./spuriousSSM score=automatic template=examples/TAE.st wc=examples/TAE.wc eq=examples/TAE.eq quiet=ALL
  ./spuriousSSM score=automatic template=examples/Oscillator.st wc=examples/Oscillator.wc eq=examples/Oscillator.eq quiet=ALL

*/


#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>


#include <math.h>
#include <ctype.h>

extern double log2(double);

#define max(a,b) ((a)>(b)?(a):(b))

#define MIN(a,b) ((a)<(b)?(a):(b))
#define MAX(a,b) ((a)>(b)?(a):(b))

// globals:

int debug=0;  int quiet=0, watch=1;
double W_verboten=0.0, W_spurious=1.0, W_bonds=0.0; 
int W_verboten_set=0, W_spurious_set=0, W_bonds_set=0;
int spurious_intraS_set=0, spurious_interS_set=0, spurious_interC_set=0, bmax_set=0;
double spurious_beta=5.0, spurious_intraS=25.0/3125, spurious_interS=5.0/3125, spurious_interC=1.0/3125, spurious_mismatch=25.0;
int spurious_equality=1, spurious_range=6;
double verboten_weak=1.0, verboten_strong=2.0, verboten_regular=0.5;
int N=0,NrS=0,NSt=0,Nwc=0,Neq=0; // note that many functions assume that their input arrays are length N
int *testwc, *testeq; char *testS, *testSt; 
char *rS_filename=NULL, *St_filename=NULL, *wc_filename=NULL, *eq_filename=NULL, *S_filename=NULL;
double temperature=37;  // default value unless overwritten by command-line option
int Nfree=0;  // number of "independent" base pairs that can be mutated, i.e. i+1 == min (eq[i],wc[i]) and St[i] isn't A, C, G, or T.
int *freeloc; // freeloc[m] is position of m^th such independent base pair (starting at 0)


// WARNING:  TREACHEROUS PROGRAMMING CONVENTION:
// wc[i], eq[i], S[i], St[i] are all indexed using 0...(N-1)
// but the values of wc[i] and eq[i] are given using 1...N, as is the case in the input files.

/*------------------------------------------------------------------------*/
/* borrowing some random number code from ViennaRNA 1.4 */
unsigned short xsubi[3];

void init_rand(void)
{
  time_t t;
  (void) time(&t);
  xsubi[0] = (unsigned short) t;
  xsubi[1] = (unsigned short) ((unsigned)t >> 16);
  xsubi[2] = 5246;
}

extern double erand48(unsigned short[]); 
double urn(void)    
     /* uniform random number generator; urn() is in [0,1] */
     /* uses a linear congruential library routine */ 
     /* 48 bit arithmetic */
{ 
  return erand48(xsubi);
}

int int_urn(int from, int to)
{
  return ( ( (int) (urn()*(to-from+1)) ) + from );
}
/*------------------------------------------------------------------------*/


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


// writes Sorig into S, as a sequence of 1,2,3,4 for A,C,G,T, or 0 if none-of-the-above
void S1234(char *S, char *Sorig)
{
  int i;  
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


// --------------------------------------------------------------------------
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
 int i, imax = (kmax-kmin+1); 
 int k,j;

 // this code could be improved by making storage static, rather than malloc'ing it at every invocation

 C = calloc(3*imax, sizeof(int)); if (C==NULL) return NULL;
 for (i=0; i<3*imax; i++) C[i] = 0;

 S   = malloc(N); if (S  ==NULL) { free(C); return NULL; }
 Swc = malloc(N); if (Swc==NULL) { free(C); free(S); return NULL; }

 S1234(S,Sorig);

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
  m = 0;              //% the ID of current subsequence
  mWC = 0;            //% the ID for this subsequence's WC complement
  for (i=1; i<=N; i++) {     //  i indexes character position in the sequence string, starting at 1
   if ( S[i-1] == 0 ) {      // reset strand & complex indicators
    Ss = i; m=0; mWC = 0;
    if ( i>1 && S[i-1 -1] == 0 ) Cs = i; 
   } else {
    m   = ( 4*m + S[i-1]-1 ) % M;  
    mWC = (mWC/4) + (Swc[i-1]-1)*M/4; 
    if (i-Ss >= k && mismatch==0) {
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
    }
    if (i-Ss >= k && mismatch==1) { int nt, bp, misWC;
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
    }
   }
  }
  C[2*imax+k-kmin+1 -1] = C[2*imax+k-kmin+1 -1] - C[imax+k-kmin+1 -1];
  C[  imax+k-kmin+1 -1] = C[  imax+k-kmin+1 -1] - C[     k-kmin+1 -1];
  free(p); free(P);
 }
 free(S); free(Swc);
 return C;
}


// --------------------------------------------------------------------------

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

// modifies S.   S MUST BE UPPERCASE ACGT & degenerate bases
void WCstring(char *S)
{
  int i,n; n=strlen(S);
  for (i=0; i<n; i++) { S[i]=WC(S[i]); }
}

// --------------------------------------------------------------------------

char degenerates[96] = "AA CC GG TT RA RG YC YT WA WT SC SG MA MC KG KT BC BG BT DA DG DT HA HC HT VA VC VG NA NC NG NT";

int test_consistency(char *S, char *St, int *wc, int *eq)
{
  int i, OK=1;
  // load_input_files guarantees that all arrays are length N

  // test conditions that must hold if wc and eq are correctly defined
  for (i=0; i<N; i++) {
    if (wc[i] != -1 && wc[wc[i]-1] != eq[i]) {
      fprintf(stderr,"ERROR: wc[wc[%d]] != eq[%d]\n",i+1,i+1); OK=0;
    }
    if (eq[i] != 0 && eq[eq[i]-1] != eq[i]) {
      fprintf(stderr,"ERROR: eq[eq[%d]] != eq[%d]\n",i+1,i+1); OK=0;
    }
  }
  if (!OK) return 0;

  // now, check that S satisfies St and wc and eq.
  for (i=0; i<N; i++) {
    char pair[3]; 
    pair[0]=St[i]; pair[1]=S[i]; pair[2]=0; // template:sequence pair for testing compatibility
    if (S[i] != ' ' && strstr(degenerates,pair)==NULL) { 
      fprintf(stderr,"ERROR: sequence S[%d]=%c is not compatible with template St[%d]=%c \n",i+1,S[i],i+1,St[i]); OK=0;
    }
    if (wc[i] != -1 && S[i] != WC(S[wc[i]-1])) {
      fprintf(stderr,"ERROR: sequence S[%d]=%c should be WC to S[%d]=%c \n",i+1,S[i],wc[i],S[wc[i]]); OK=0;
    }
    if (eq[i] != 0 && S[i] != S[eq[i]-1]) {
      fprintf(stderr,"ERROR: sequence S[%d]=%c should be equal to S[%d]=%c \n",i+1,S[i],eq[i],S[eq[i]]); OK=0;
    }
  }
  // now, check that St satisfies wc and eq.
  for (i=0; i<N; i++) {
    if (wc[i] != -1 && St[i] != WC(St[wc[i]-1])) {
      fprintf(stderr,"ERROR: template St[%d]=%c should be WC to St[%d]=%c \n",i+1,St[i],wc[i],St[wc[i]]); OK=0;
    }
    if (eq[i] != 0 && St[i] != St[eq[i]-1]) {
      fprintf(stderr,"ERROR: template St[%d]=%c should be equal to St[%d]=%c \n",i+1,St[i],eq[i],St[eq[i]]); OK=0;
    }
  }
  return OK;
}

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
  int i,j; int *marked;

  marked = (int *) calloc(N, sizeof(int));

  for (i=0; i<N; i++) {
    if (!marked[i]) {
      // class = find(eq==eq(i));
      // for j=1:length(class), S(class(j))=S(i); end
      // marked(class)=1;

      for (j=0; j<N; j++) {
        if (eq[j] == eq[i]) { S[j] = S[i]; marked[j]=1; }
      }

      // class = find(eq==wc(i));
      // for j=1:length(class), S(class(j))=WC(S(i)); end
      // marked(class)=1;

      for (j=0; j<N; j++) {
        if (eq[j] == wc[i]) { S[j] = WC(S[i]); marked[j]=1; }
      }
    }
  }
  free(marked);

}

void constrain_single_fast(char *S, int *wc, int *eq, int i)
{
  int j;

  for (j=i+1; j<N; j++) {
    if (eq[j] == eq[i]) { S[j] = S[i];}
  }

  for (j=i+1; j<N; j++) {
    if (eq[j] == wc[i]) { S[j] = WC(S[i]); }
  }

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
%   DNAdesign-0.02          Erik Winfree
*/

//  This was the old DNAdesign and spuriousC scoring.  See help info for revised spuriousSSm scoring that is actually implemented below.
//    "   the sum over all 6-nucleotide windows of:\n"
//    "      how many of these patterns match the window: WWWWWW, SSSSSS, RRRRRR, YYYYYY, RYRYRY, YRYRYRYR, SWWWWW, WWWWWS, SSSSSW, WSSSSS\n"
//    "      + how many substrings appear in the window: GGG, GGGG, CCC, CCCC, TTTT, AAAA, GCGC, GGCC, CCGG, CGCG\n"
//    "      + 0.1 * (there are 4 or more W bases in the window).\n"


double bad6mer[1<<12];

void make_bad6mer()   
{
  char S[7]; char ACGT[5]=" ACGT";
  int i,j, b1,b2,b3,b4,b5,b6;
  int co[5],ce[5],c[5];

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
         bad6mer[i] =  // these penalties have evolved since spuriousC...
	   verboten_weak * (
			    (c[2]+c[3]==0) +             // all AT
			    (strstr(S,"TTTT")!=NULL) + 
			    (strstr(S,"AAAA")!=NULL) + 
			    ((c[1]+c[4])==5 && (b1==2 || b1==3)) +  // 5 AT is weak
			    ((c[1]+c[4])==5 && (b6==2 || b6==3)) +  // 5 AT is weak
			    (c[1]+c[4]>3)                           // 4 or more AT within 6 is weak
	   ) + verboten_strong * (
			    (c[1]+c[4]==0) +             // all CG
			    (strstr(S,"GGG")!=NULL) + 
			    1000*(strstr(S,"GGGG")!=NULL) + 
			    (strstr(S,"CCC")!=NULL) + 
			    1000*(strstr(S,"CCCC")!=NULL) + 
			    ((c[2]+c[3])==5 && (b1==1 || b1==4)) +  // 5 CG is bad
			    ((c[2]+c[3])==5 && (b6==1 || b6==4))    // 5 CG is bad
	   ) + verboten_regular * (
			    (c[1]+c[3]==0) +                // all AG (Pu)
			    (c[2]+c[4]==0) +                // all TC (Py)
			    (co[1]+co[3]+ce[2]+ce[4]==0) +  // alt PuPy
			    (ce[1]+ce[3]+co[2]+co[4]==0) +  // alt PyPu
			    (strstr(S,"GCGC")!=NULL) + 
			    (strstr(S,"GGCC")!=NULL) + 
			    (strstr(S,"CCGG")!=NULL) + 
			    (strstr(S,"CGCG")!=NULL)  
	   );
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
%   DNAdesign-0.02          Erik Winfree
*/

// note: could merge with S1234() defined above, which converts entire array
static inline int C1234(char c)   //  map A,C,G,T to 1,2,3,4
{ 
  int v=0;
  switch (c) {
    case 'A': case 'a': case 1: v=1; break;
    case 'C': case 'c': case 2: v=2; break;
    case 'G': case 'g': case 3: v=3; break;
    case 'T': case 't': case 4: v=4; break;
    default: v=0; 
  }
  return v;
}

double score_verboten(char *S, int *wc, int *eq)
{
 int i, K, M, m; double V;

 V=0.0; K=5;

 M = 1<<12;
 m = 0;                //% the ID of current subsequence
 for (i=0; i<N; i++) {
  if (C1234(S[i]) == 0) {     //% reset strand & complex indicators
    m=0; K=5;          //% wait till full 6-mer considered
  } else {
    m = (  (m<<2) + C1234(S[i])-1 ) & (M-1);
    if (K) K=K-1; else V=V+bad6mer[m]; 
  }				       
 }

 // printf("verboten('%s')=%f\n\n",S,V);

 return V; 
}


// optimizing score_spurious performs so poorly on Milo's test;
// Niles & Robert's positive+negative approach is more successful...
// add some element of positive design...
//    here, we sum the nearest-neighbor Hbond/stacking energy for
//    all desired WC interactions, in Kcal/mole. 
double score_bonds(char *S, int *wc, int *eq)
{
  double score=0; int i; 

  //  for (i=0; i<5; i++) 
  //     { for (N=0; N<5; N++) printf("%4.2f ",nnDH[i][N]); printf("\n"); }

  for (i=1; i<N; i++)
    if (wc[i-1]==wc[i]+1) score += nnDG(C1234(S[i-1]),C1234(S[i]));

  return score/2;  // every nt involved in WC get counted twice per stack 
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
  switch(Stc) { 
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
  default:  choices = " ";
  }
  return choices[ int_urn(0,strlen(choices)-1) ];
}

// make a single random change, apply the constrains, and insist that the change wasn't eliminated.
// doesn't make sense to call this if Nfree==0, in which case no mutation is made.
void mutate(char *S, char *St, int *wc, int *eq) 
{
  int i; char oldc;
  if (Nfree==0) return;
  i = freeloc[int_urn(0,Nfree-1)];
  oldc = S[i];
  do {
    S[i] = randbasec(St[i]);   
  } while (oldc == S[i]);
  constrain_single_fast(S,wc,eq,i);
}
// --------------------------------------------------------------------------

int set_auto_spurious_weights(char *S, int *wc, int *eq) // returns bmax value
{
  // n_i are strand lengths;  c are complexes, sets of strands
  // weight for length-k intra-strand bindings is beta^(k-k0) where k0 = log_4 (max_i n_i)  
  // weight for length-k inter-strand intra-complex bindings is beta^(k-k1) where k1 = log_4 (max_c sum_{i in c} n_i)
  // weight for length-k inter-complex bindings is beta^(k-k2) where k2 = log_4 (sum_i n_i)
  // the idea here being that k0, k1, k2 represent the total number of unique subsequences needed

  int i, n, nq=0, nk0=0, nk1=0, nk2=0, ni=0, nc=0, nbp=0;

  i=0; while (S[i]!=0) {
    if (S[i] != ' ') { nk2++; ni++; nc++; nk0=max(nk0,ni); nk1=max(nk1,nc); }
    if (S[i] == ' ') { ni=0; }
    if (S[i] == ' ' && i>0 && S[i-1] == ' ') { nc=0; }
    i++;
  }

  if (!spurious_intraS_set) spurious_intraS = pow(spurious_beta,0.0-log2(1.0*nk0)/2);
  if (!spurious_interS_set) spurious_interS = pow(spurious_beta,0.0-log2(1.0*nk1)/2);
  if (!spurious_interC_set) spurious_interC = pow(spurious_beta,0.0-log2(1.0*nk2)/2);

  if (!W_spurious_set) W_spurious = 1.0;            // logic: expect a score of about 1.0 for each subcategory
  if (!W_verboten_set) W_verboten = 640.0 / nk2;     // logic: three-letter verboten sequences occur once/64, so "typical" score should be 10.
  
  for (i=1; i<N; i++) if (wc[i-1]==wc[i]+1) nbp++;  // count number of defined base-pairs 
  if (!quiet) printf("Automatic: counted %d base-pairing stacks in target structures.\n",nbp);
  if (nbp==0) nbp++;
  if (!W_bonds_set) W_bonds = 10.0/nbp;             // logic:  each stack varies from about 1.0 to 2.0 kcal/mol, so let's go for "typical" score 10.

  // count number of unique base equivalent classes
  n=strlen(S); for (i=0; i<n; i++)  if (eq[i]==i+1 && (wc[i]>i+1 || wc[i]==-1)) nq++; 
  if (!quiet) printf("Automatic: counted %d unique base equivalence classes.\n",nq);

  return (12*nq+1); // get bored after 4-fold coverage of all possible single mutations
}

// scoring function from score_spurious.m in DNAdesign
double score_spurious(char *S, int *wc, int *eq)
{
  // prefactors = [5^(-4) 5^(-5) 5^(-6)];  % most weight to intramolecular
  // S0 = prefactors * (spurious(S, 3,8, wc,eq) * cumprod(5*ones(1,6))')  ;
  // S1 = prefactors * (spurious1(S, 5,10) * cumprod(5*ones(1,6))')  ;
  
  double *factors;
  int *C; double score=0.0; int i, imax=3*spurious_range;

  factors = (double *) malloc(imax*sizeof(double));  // typically spurious_range=6 and there are 18 factors.
  for (i=0; i<imax; i++) factors[i]=0.0;

  factors[0*spurious_range] = spurious_intraS; for (i=1; i<spurious_range; i++) factors[i]                  = factors[i-1]                  *spurious_beta;
  factors[1*spurious_range] = spurious_interS; for (i=1; i<spurious_range; i++) factors[i+1*spurious_range] = factors[i-1+1*spurious_range] *spurious_beta;
  factors[2*spurious_range] = spurious_interC; for (i=1; i<spurious_range; i++) factors[i+2*spurious_range] = factors[i-1+2*spurious_range] *spurious_beta;

  C = spurious(0,S,3,2+spurious_range,wc,eq);
  for (i=0; i<imax; i++) score+=factors[i]*C[i];
  free(C);

  C = spurious(1,S,5,4+spurious_range,NULL,NULL);
  for (i=0; i<imax; i++) score+=factors[i]*C[i]*spurious_mismatch/spurious_beta/spurious_beta;
  free(C);

  free(factors);
  return score;
}
// --------------------------------------------------------------------------


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



double score_all(char *rS, int *wc, int *eq)
{
  double V=0.0, Spu=0.0, Bnd=0.0;

  if (W_verboten>0) V   = score_verboten(rS,wc,eq);
  if (W_spurious>0) Spu = score_spurious(rS,wc,eq);
  if (W_bonds>0)    Bnd = score_bonds(rS,wc,eq);

  return W_verboten*V + W_spurious*Spu + W_bonds*Bnd; 
}


// --------------------------------------------------------------------


void test_evals(char *testS, int *testwc, int *testeq)
{
 int *C; int i,k; double V=0, Spu=0, Bnd=0;
 int sh = (quiet==0 && watch==0);

 printf("\n");

 if (W_spurious>0) {

   C = spurious(0, testS, 3,2+spurious_range, testwc, testeq);
   if (spurious_equality) printf("spurious counts identity matches as well as WC matches.\n");
   else printf("spurious counts WC matches but not identity matches.\n");
   printf("spurious(testS, 3,%d, testwc, testeq)\n",2+spurious_range);
   printf("C =  \n");
   for (i=0;i<3;i++) {
     for (k=3; k<=2+spurious_range; k++) {
       printf(" %5d", C[i*spurious_range + k-3]);
     }
     printf("\n");
   }
   free(C);

   C = spurious(0, testS, 3,2+spurious_range, NULL, NULL);
   printf("spurious0(testS, 3,%d)\n",2+spurious_range);
   printf("C0 =  \n");
   for (i=0;i<3;i++) {
     for (k=3; k<=2+spurious_range; k++) {
       printf(" %5d", C[i*spurious_range + k-3]);
     }
     printf("\n");
   }
   free(C);


   C = spurious(1, testS, 5,4+spurious_range, NULL, NULL);
   printf("spurious1(testS, 5,%d)\n",4+spurious_range);
   printf("C1 =  \n");
   for (i=0;i<3;i++) {
     for (k=5; k<=4+spurious_range; k++) {
       printf(" %5d", C[i*spurious_range + k-5]);
     }
     printf("\n");
   }
   free(C);
   printf("spurious: intraS = %11.5f, interS = %11.5f, interC  = %11.5f, beta = %5.3f, mismatch = %5.3f\n", 
	  spurious_intraS, spurious_interS, spurious_interC, spurious_beta, spurious_mismatch);
 }
 if (W_verboten>0) printf("verboten: weak   = %11.5f, strong = %11.5f, regular = %11.5f\n", verboten_weak, verboten_strong, verboten_regular);

 if (sh || W_verboten>0) printf("** score_verboten score = %11.5f \n", V=score_verboten(testS,testwc,testeq) );
 if (sh || W_spurious>0) printf("** score_spurious score = %11.5f \n", Spu=score_spurious(testS,testwc,testeq) );
 if (sh || W_bonds>0) printf("** score_bonds    score = %11.5f \n", Bnd=score_bonds(testS,testwc,testeq) );
 printf("** [verboten spurious bonds] = [%11.5f %11.5f %11.5f]-weighted score = %11.5f \n", 
	W_verboten, W_spurious, W_bonds, 
	W_verboten*V + W_spurious*Spu + W_bonds*Bnd);

 printf("\n\n");

}



// -------------------------------------------------------------------------

long int get_file_length(char *name)
{
  FILE *fp;
  int len;

  if (name==NULL) return 0;

  fp = fopen(name, "r");
  if( fp == NULL ) {
    fprintf(stderr,"Error opening file %s.", name); exit(-1);
  }
  fseek(fp, 0, SEEK_END);
  len = ftell(fp);
  fclose(fp);
  return len;
}

void load_input_files()
{ FILE *f; char c; int i; double r; int maxN;

//// Note: we strip off trailing ' ' from rS and St,
//// and delete corresponding entries in wc and eq
//// but YELL if any deleted entries are not -1 or 0 respectively.

 maxN = 100+get_file_length(rS_filename);
 testS = malloc(maxN);
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

 maxN = 100+get_file_length(St_filename);
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

 maxN = 100+get_file_length(wc_filename);  
 testwc = calloc(maxN, sizeof(int));  // plenty of space, since ints take more than 1 char
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

 maxN = 100+get_file_length(eq_filename);
 testeq = calloc(maxN, sizeof(int)); // plenty of space, since ints take more than 1 char
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
   free(testSt); testSt = malloc(N+100);
   NSt=N; for (i=0; i<N; i++) testSt[i]='N'; testSt[N]=0;
   if (debug) { fprintf(stderr,"St = <"); for (i=0; i<NSt; i++) fprintf(stderr,"%c",testSt[i]); fprintf(stderr,">\n"); }
 }
 if (rS_filename==NULL) {
   free(testS); testS = malloc(N+100); 
   NrS=NSt; for (i=0; i<NrS; i++) testS[i]=randbasec(testSt[i]); testS[NrS]=0;
   if (debug) { fprintf(stderr,"S = <"); for (i=0; i<NrS; i++) fprintf(stderr,"%c",testS[i]); fprintf(stderr,">\n"); }
 }
 if (wc_filename==NULL) {
   free(testwc); testwc = calloc(N+100, sizeof(int)); 
   Nwc=N; for (i=0; i<N; i++) testwc[i]=-1;
   if (debug) { fprintf(stderr,"wc = <"); for (i=0; i<Nwc; i++) fprintf(stderr,"%d ",testwc[i]); fprintf(stderr,">\n"); }
 }
 if (eq_filename==NULL) {
   free(testeq); testeq = calloc(N+100, sizeof(int)); 
   Neq=N; for (i=0; i<N; i++) testeq[i]=(testS[i]==' ')?0:(i+1); 
   if (debug) { fprintf(stderr,"eq = <"); for (i=0; i<Neq; i++) fprintf(stderr,"%d ",testeq[i]); fprintf(stderr,">\n"); }
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

}

// -------------------------------------------------------------------------
int main(int argc, char *argv[])
{
  int i; FILE *f;
  time_t t_start,t_now,tmax=0;  int imax=0; int bored, bmax=0; int automatic=0;
  int fancy_args=1, trace=0; 

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
    "  sequence=[file]      for initial nucleotide sequence (rS_file) [default: random]\n"
    "  template=[file]      for allowed nucleotides at each position (St_file) [default: all N]\n"
    "  wc=[file]            for Watson-Crick pairing constraints (wc_file) [default: all -1]\n"
    "  eq=[file]            for equality constraints (eq_file) [default: identity]\n"
    "  output=[filename]    for final output sequence upon completion [default: stdout]\n"
    "  N=[value]            length of sequence, if no files are specified\n"
    "  quiet=TRUE           don't print anything, except errors & output\n"
    "  quiet=SCORES         output intial & final test scores, for all relevant score methods\n"
    "  quiet=WATCH          output running tally of best test scores, for chosen method\n"
    "  quiet=ALL            output all the above [default]\n"
    "  trace=ON             output the current sequence, at every improvement\n"
    " STOPPING CRITERIA\n"
    "  tmax=[value]         stop after specified time (in seconds) has elapsed\n"
    "  imax=[value]         stop after specified number of mutation iterations have been completed\n"
    "  bmax=[value]         stop after specified number of mutation iterations yield no change\n"
    " SCORING FUNCTION\n"
    "  score=automatic      set stopping criteria and spurious, verboten, and bonds score weights based on design size\n"
    "  score=noautomatic    reset any previous score=automatic option that has been provided.\n"
    "  score=[verboten|spurious|bonds] what scoring function to use [default: spurious]\n"
    "     verboten: penalizes use of subsquences such as GGG, AAAA, WWWWWW, PuPyPuPyPuPy\n"
    "     spurious: penalizes pairs of undesired exact WC matches of lengths 3-8 (and 1-mm of lengths 5-10)\n"
    "     bonds:    prefers strong nearest-neighbor stacking energies for designed helices\n"
    "     setting 'score=xxx' is shorthand for 'W_xxx=1.0' and other weights 0.0; this can be further adjusted...\n"
    "  W_verboten=[value]               weight for verboten score [default: 0.0]\n"
    "  W_spurious=[value]               weight for spurious score [default: 1.0]\n"
    "  W_bonds=[value]                  weight for bonds score [default: 0.0]\n"
    "  spurious_beta=[value]            multiplicative penalty, per match base pair [default: 5.0]\n"
    "  spurious_intraS=[value]          initial penalty, intra-strand matches [default: 1/125]\n"
    "  spurious_interS=[value]          initial penalty, inter-strand matches [default: 1/625]\n"
    "  spurious_interC=[value]          initial penalty, inter-complex matches [default: 1/3125]\n"
    "  spurious_mismatch=[value]        mismatch penalty [default: 25]\n"
    "  spurious_range=[value]           number of values of k (match length) tested [default: 6]\n" 
    "  spurious_equality=[0|1]          penalize identity matches, not just WC complementary matches [default: 1]\n"
    "  verboten_weak=[value]            multiplicative penalty for weak subsequences [default: 1.0]\n"
    "  verboten_strong=[value]          multiplicative penalty for strong subsequences [default: 2.0]\n"
    "  verboten_regular=[value]         multiplicative penalty for regular subsequences [default: 0.5]\n"
    "  temperature=[value]              temperature in Celsius for bond energy calculations [default: 37]\n"
    "\n"
    "Notes:\n\n"
    "sequence is in { A C G T }*\n\n"
    "template is in { A C G T  R  Y  W  S  M  K  B   D   H   V   N   }*\n"
    "        meaning  A C G T AG CT AT CG AC GT CGT AGT ACT ACG ACGT \n\n"
    "        spaces indicate strand (' ') and complex ('  ') boundaries\n\n"
    "The eq() and wc() arrays are lists of integers, indexes into the sequence string, starting at 1. \n"
    "   eq(i) is the min of all bases constrained to have the same value as i (' ' bases have eq(i)==0).\n"
    "   wc(i) is the min of all bases constrained to be complementary to i    (or -1 if there are none).\n"
    "Strand id's are numbered starting at 1.\n\n"
    "It is the responsibility of the user to ensure that indirect implications of wc() and eq() are considered, \n"
    "   e.g. the complement of a complement should be considered equal, and that equality is transitive.\n"
    "   If wc() and eq() do not uphold these properties, odd things might happen.  Beware.\n"
    "\n"
    "The spurious score:\n"
    "In the table C = spurious(S, kmin, kmax, wc, eq)\n"
    "   where 1 <= i <= 3 (intramolec, intracomplex, intercomplex match) \n"
    "   and kmin <= k <= kmax (match length), \n"
    "C(i,k) gives the number of undesired WC subsequences of length k \n"
    "   that match another subsequence exactly or are exactly complementary [modulated by spurious_equality], \n"
    "   where we exclude matches that are fully implied by the wc() and eq() inputs \n"
    "   and we count self-complementary as well as self-overlapping matches; \n"
    "   matches are divided according to how ''local'' the match is: \n"
    "       i=1  intramolecular              (same strand) matches \n"
    "       i=2  intermolecular intracomplex (different strand, same complex) matches \n"
    "       i=3  intermolecular intercomplex (different strand, different complex) matches \n"
    "The related table C0 = spurious0(S, kmin, kmax)\n"
    "   provides the same information, but doesn't exclude matches that are implied by wc() and eq().\n"
    "   This table is not used in scoring, but printed as information potentially useful to the user.\n"
    "The related table C1 = spurious1(S, kmin, kmax)\n"
    "   provides similar information, but only counts matches that have exactly 1 mismatch.\n"
    "   Bulges (e.g. resulting from insertions or deletions) are not considered.\n"
    "The spurious score is:\n"
    "   sum_{i in 1=intraS...3=interC},{k in 3..8}  C(i,k) * beta^(k-3) * spurious_i\n"
    "   + sum_{i in 1=intraS...3=interC},{k in 5..10}  C1(i,k) * beta^(k-5) * spurious_i * spurious_mismatch\n"
    "The upper value for k will increase or decrease if spurious_range is changed from its default value of 6.\n"
    "\n"
    "The verboten score is:\n"
    "   the sum over all 6-nucleotide windows, of penalties for whether the window\n"
    "       exactly matches certain patterns, contains certain substrings, or has certain base counts: \n"
    "   verboten_weak    * ( #match of {WWWWWW, SWWWWW, WWWWWS} + #contained of {TTTT, AAAA} + #if {has 4 or more W bases} ) + \n"
    "   verboten_strong  * ( #match of {SSSSSS, SSSSSW, WSSSSS} + #contained of {GGG, GGGG, CCC, CCCC} ) + \n"
    "   verboten_regular * ( #match of {RRRRRR, YYYYYY, RYRYRY, YRYRYR} + #contained of {GCGC, GGCC, CCGG, CGCG} ) \n"
    "       where the GGGG and CCCC matches are multiplied by 1000 to ensure those sequences are eliminated. \n"
    "\n"
    "The bonds score is:\n"
    "   the sum over all base-pair stacks, according to wc(), of the SantaLucia 1998 nearest-neighbor free energy.\n"
    "   Uses the temperature parameter and dG = dH - T dS; stronger bonding gives more negative values, i.e. improves the score.\n"
    "\n"
    "The overall score is W_spurious * score_spurious + W_verboten * score_verboten + W_bonds * score_bonds.\n\n"
    "   Thus, weighted combinations of scoring functions can be specified using the 'W_*' parameters.\n");

    exit(-1);
 }

 init_rand();

 if (fancy_args) {
   // try to parse fancy args

   for (i=1;i<argc;i++) {
     if      (strncmp(argv[i],"sequence=",9)==0) rS_filename=&argv[i][9]; 
     else if (strncmp(argv[i],"template=",9)==0) St_filename=&argv[i][9]; 
     else if (strncmp(argv[i],"wc=",3)==0) wc_filename=&argv[i][3]; 
     else if (strncmp(argv[i],"eq=",3)==0) eq_filename=&argv[i][3]; 
     else if (strncmp(argv[i],"output=",7)==0) S_filename=&argv[i][7]; 
     else if (strncmp(argv[i],"N=",2)==0) N=atoi(&argv[i][2]);
     else if (strncmp(argv[i],"tmax=",5)==0) tmax=atoi(&argv[i][5]);
     else if (strncmp(argv[i],"imax=",5)==0) imax=atoi(&argv[i][5]);
     else if (strncmp(argv[i],"bored=",6)==0) {bmax=atoi(&argv[i][6]); bmax_set=1;}
     else if (strncmp(argv[i],"bmax=",5)==0) {bmax=atoi(&argv[i][5]); bmax_set=1;}
     else if (strncmp(argv[i],"score=automatic",15)==0) automatic=1;
     else if (strncmp(argv[i],"score=noautomatic",17)==0) automatic=0;
     else if (strncmp(argv[i],"score=verboten",14)==0) { W_spurious=0.0; W_verboten=1.0; W_bonds=0.0; }
     else if (strncmp(argv[i],"score=spurious",14)==0) { W_spurious=1.0; W_verboten=0.0; W_bonds=0.0; }
     else if (strncmp(argv[i],"score=bonds",11)==0) { W_spurious=0.0; W_verboten=0.0; W_bonds=1.0; }
     else if (strncmp(argv[i],"quiet=TRUE",10)==0) { quiet=1; watch=0; }
     else if (strncmp(argv[i],"quiet=SCORES",12)==0) { quiet=0; watch=0; }
     else if (strncmp(argv[i],"quiet=WATCH",11)==0) { quiet=1; watch=1; }
     else if (strncmp(argv[i],"quiet=ALL",9)==0) { quiet=0; watch=1; }
     else if (strncmp(argv[i],"trace=ON",8)==0) { trace=1;  }
     else if (strncmp(argv[i],"debug=TRUE",10)==0) debug=1;
     else if (strncmp(argv[i],"W_verboten=",11)==0) {W_verboten=atof(&argv[i][11]); W_verboten_set=1;}
     else if (strncmp(argv[i],"W_spurious=",11)==0) {W_spurious=atof(&argv[i][11]); W_spurious_set=1;}
     else if (strncmp(argv[i],"W_bonds=",8)==0) {W_bonds=atof(&argv[i][8]); W_bonds_set=1;}
     else if (strncmp(argv[i],"spurious_beta=",14)==0) spurious_beta=atof(&argv[i][14]);
     else if (strncmp(argv[i],"spurious_intraS=",16)==0) {spurious_intraS=atof(&argv[i][16]); spurious_intraS_set=1;}
     else if (strncmp(argv[i],"spurious_interS=",16)==0) {spurious_interS=atof(&argv[i][16]); spurious_interS_set=1;}
     else if (strncmp(argv[i],"spurious_interC=",16)==0) {spurious_interC=atof(&argv[i][16]); spurious_interC_set=1;}
     else if (strncmp(argv[i],"spurious_mismatch=",18)==0) {spurious_mismatch=atof(&argv[i][18]); }
     else if (strncmp(argv[i],"spurious_range=",15)==0) {spurious_range=atoi(&argv[i][15]); }
     else if (strncmp(argv[i],"spurious_equality=",18)==0) spurious_equality=(atoi(&argv[i][18])>0);
     else if (strncmp(argv[i],"verboten_weak=",14)==0) {verboten_weak=atof(&argv[i][14]);}
     else if (strncmp(argv[i],"verboten_strong=",16)==0) {verboten_strong=atof(&argv[i][16]);}
     else if (strncmp(argv[i],"verboten_regular=",17)==0) {verboten_regular=atof(&argv[i][17]);}
     else if (strncmp(argv[i],"temperature=",12)==0) temperature=atof(&argv[i][12]);
     else { 
       fprintf(stderr,"Cannot parse argument <%s>.  Aborting. \n\n",argv[i]); exit(-1); 
     }
   }
 
 } else { 
   // otherwise, try to read all four input files
   rS_filename=argv[1]; St_filename=argv[2]; 
   wc_filename=argv[3]; eq_filename=argv[4];
 }
 if (spurious_range < 1) spurious_range=1;
 if (spurious_range > 10) spurious_range=10;

 make_bad6mer();

 // load in the files!
 load_input_files();

 //////////////////////////////////////////////////////////////
 // information was read in.  now do the designing!

 if (automatic) {
   if (!bmax_set) {
     bmax=set_auto_spurious_weights(testS,testwc,testeq);
   }
   else set_auto_spurious_weights(testS,testwc,testeq);
 }

 constrain(testS,testwc,testeq);
 if (!test_consistency(testS,testSt,testwc,testeq)) {
     fprintf(stderr,"ERROR: input files have non-explicit implications.\n\n"); exit(-1); 
 }

 if (!quiet) printf("\nconstrained S = <%s>  N=%d \n",testS,N);

 freeloc = calloc( N, sizeof(int) );
 Nfree = 0;
 for (i=0; i<N; i++) {
     if ( (testeq[i]==i+1) && (testwc[i]>i+1 || testwc[i]==-1) && (testSt[i]!='A' && testSt[i]!='C' && testSt[i]!='G' && testSt[i]!='T' ) ) { 
         freeloc[Nfree]=i; Nfree++;
     }
 }
 if (!quiet) printf("\nFound %d bases that can probably be changed freely.\n\n",Nfree);

 if (!quiet) test_evals(testS,testwc,testeq);



 {
  char *oldS; double old_score, new_score; int steps=0; bored=0;
  oldS = malloc(N); 
  old_score = score_all(testS,testwc,testeq); 
  if (watch) printf("%8d steps, %8d seconds : score = %18.10f",0,0,old_score);
  if (watch && (imax>0)) printf(" (imax=%d)", imax);
  if (watch && (tmax>0)) printf(" (tmax=%d)", (int) tmax);
  if (watch && (bmax>0)) printf(" (bored=%d,bmax=%d)", bored, bmax);
  if (watch) printf("\n");
  
  time(&t_start);  time(&t_now);
  while ( (tmax==0 || t_start+tmax>=t_now) && (imax==0 || steps<imax) && (bmax==0 || bored<bmax) && Nfree>0) { 
   steps++; 
   if (!(quiet && !watch) && !trace && steps%300 == 0 && tmax==0 && imax==0 && bmax==0) 
      printf("\n\n%s\n\n",testS);
   for (i=0; i<N; i++) oldS[i]=testS[i];
   mutate(testS,testSt,testwc,testeq);
   // if (bored > bmax/2) mutate(testS,testSt,testwc,testeq); // get desperate and start mutating two at a time (doesn't seem to work)
   new_score = score_all(testS,testwc,testeq); 
   if (new_score <= old_score) {
     if (watch) printf("%8d steps, %8d seconds : score = %18.10f", steps, (int)(t_now-t_start), new_score);
     if (watch && (imax>0)) printf(" (imax=%d)", imax);
     if (watch && (tmax>0)) printf(" (tmax=%d)", (int) tmax);
     if (watch && (bmax>0)) printf(" (bored=%d,bmax=%d)", bored, bmax);
     if (watch) printf("\n");
     if (trace) printf("%s\n",testS);
     if (new_score < old_score) bored=0;
     old_score = new_score; 
     fflush(stdout);
   } else {
     for (i=0; i<N; i++) testS[i]=oldS[i]; bored++;
   }
   time(&t_now);
  }  
  if (watch) printf("%8d steps, %8d seconds : score = %18.10f FINAL\n", 
               steps-1, (int)(t_now-t_start), old_score);
 }

 //// print final sequence to output file *****

 constrain(testS,testwc,testeq);
 if (!test_consistency(testS,testSt,testwc,testeq)) {
     fprintf(stderr,"ERROR: output has non-explicit implications! THIS SHOULDN'T HAPPEN!\n\n"); exit(-1); 
 }

 if (!quiet) test_evals(testS,testwc,testeq);


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

 return 0;

}
