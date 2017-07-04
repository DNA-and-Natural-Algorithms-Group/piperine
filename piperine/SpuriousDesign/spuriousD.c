// SpuriousD++ by David Zhang
// based off SpuriousC by Erik Winfree

// Last Modified 7/23/2001, DZ.  

// 7/23/2001: Fixed "bulge" algorithms to include template side bulges
//            as well as sense side bulges.
// 7/16/2001: Added "bounce-up" functionality.
//            Code is less likely to be trapped in local minima now.
// 7/2/2001:  Modified parameters.
//            Added G-T wobble size 2 bulges.
// 6/27/2001: Added G-T wobble bulges.
//            Seems to work well.  
// 6/26/2001: Debugged G-T non-mismatched.  Seems to work.
//            Also added bulges to spurious1.  DZ
// 6/25/2001: Added G-T as base-pair for purposes of scores, but not for 
//            restrictions via WC file.  DZ


/*  Converted to C from Matlab by hand, EW, 4/00 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

#define max_GT_wobbles 3

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
*/

// UpdateGTScores Function.  Updates all scores with G-T wobble included.  Recursive function that calls itself.
// Capped at max_gt_wobbles recursions to prevent 2^n slowdown.  

void UpdateGTScores(int* C, int* p, int* P, int mWC, int k, int kmin, int imax, int Ss, int Cs, int recur, int pos, int opt)
{
  int iGT, j;
  iGT = mWC>>(2 * pos);
  //if (recur == 0)
    //    printf("\nHit %d %d %d %d",mWC, iGT, pos, recur);
  if (iGT == 0) {
    if ((recur != 0)&&(opt == 0)) {       // there was at least one G-T pairing
      j = P[mWC];
      while (j != 0) {
	if (j > Ss) C[k-kmin]++;
	if (j > Cs) C[imax+k-kmin]++;
	C[2*imax+k-kmin]++;
	j = p[j-1];
      }
    }
    else if (opt == 1) {
      j = P[mWC];
      while (j != 0) {
	if (j > Ss) C[k-kmin]++;
	if (j > Cs) C[imax+k-kmin]++;
	C[2*imax+k-kmin]++;
	j = p[j-1];
      }
    }
  }
  else {
    if (((iGT % 4 == 0)||(iGT % 4 == 1)) && (recur < max_GT_wobbles))                     // limit to max_GT_wobbles per subsequence
      UpdateGTScores(C, p, P, ((2<<(2*pos)) + mWC), k, kmin, imax, Ss, Cs, recur+1, pos+1, opt);
    UpdateGTScores(C, p, P, mWC, k, kmin, imax, Ss, Cs, recur, pos+1, opt);
  }
}


int *spurious(int mismatch, char *Sorig, int kmin, int kmax, int *wc, int *eq)
{
 int *C; char *S,*Swc;
 int i, imax = (kmax-kmin+1); 
 int N; int k,j;

 C = calloc(3*imax, sizeof(int)); if (C==NULL) return NULL;
 for (i=0; i<3*imax; i++) C[i] = 0;

 N = strlen(Sorig);  // note: wc & eq must be length N, or ELSE!!! 

 S   = malloc(N); if (S  ==NULL) { free(C); return NULL; }
 Swc = malloc(N); if (Swc==NULL) { free(C); free(S); return NULL; }

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
 for (i=0; i<N; i++) {                   //% WC bases, same order
   if (S[i]==0) Swc[i]=0; else Swc[i]=5-S[i]; 
 }

 for (k=kmin; k<=kmax; k++) { 
  int M, *P, *p, Cs, Ss, m, mWC, tempvar;
  M = 1<<(2*k);  //% 4^k;
  P = calloc(M,sizeof(int));  //% hold index of present subsequences
  p = calloc(N,sizeof(int));  //% holds pointers to more instances of sequence
  Cs = 0;             //% the start of most recent complex
  Ss = 0;             //% the start of most recent strand
  m = 0;              //% the encoding of current subsequence in binary
  mWC = 0;            //% the encoding for this subsequence's WC complement in binary
  for (i=1; i<=N; i++) {
   if ( S[i -1] == 0 ) {      //% reset strand & complex indicators
    Ss = i; m=0; mWC = 0;
    if ( i>1 && S[i-1 -1] == 0 ) Cs = i; 
   } else {
    m   = ( 4*m + S[i -1]-1 ) % M;        // sub-portion of DNA that we are looking at
    mWC = (mWC/4) + (Swc[i -1]-1)*M/4;    // Exact complement to subportion m.
    if (i-Ss >= k && mismatch==0) {
       p[i -1]=P[m]; P[m]=i; //% add site to list of subsequences

       UpdateGTScores(C, p, P, mWC, k, kmin, imax, Ss, Cs, 0, 0, 0);    // Count the matches with G-T wobbles
	 
	   
       j=P[mWC];     //% look for previous occurrences, including self
       while (j != 0) { int msum,jjj;                                   // Count matches without G-T wobbles
         if (wc != NULL && eq != NULL) { //% exclude expected match
           for (msum = 0, jjj=0; jjj<k; jjj++)
              if (wc[i-k+1 +jjj-1] == eq[j -jjj-1]) msum++;
	 } else msum=0;
         if (msum != k) {                      
           if (j>Ss) C[       k-kmin+1 -1]++;
           if (j>Cs) C[imax  +k-kmin+1 -1]++;
	             C[2*imax+k-kmin+1 -1]++;
         }
         j = p[j -1];
       }
    }

    // one-mismatches (one-bulge on each side), including G-T wobbles

    if (i-Ss >= k && mismatch==1) { int nt, bp, bp2, misWC;
     for (nt=2; nt<=k-1; nt++) {  //% where is mismatch
      misWC = mWC & (-1 - (3<<(2*nt-2)));
      for (bp=0; bp<4; bp++) { //% what is value at mismatch
       misWC += bp<<(2*nt-2);


       if (misWC != mWC) {         // has to actually be a mismatch (to prevent double counting)
	 tempvar = (mWC>>(2*nt-2)) % 4;     // initial value of replaced base
	 // exclude perfect match with G-T wobble:
	 if ((tempvar == 3)||(tempvar == 2)||((tempvar == 1) && (bp != 3))||((tempvar == 0)&&(bp != 2)))
	   UpdateGTScores(C, p, P, misWC, k, kmin, imax, Ss, Cs, 0, 0, 1);  // updates scores of one-mismatches
       }
       misWC -= bp<<(2*nt-2);
      }
     }
    
     // one-bulges, sense-side, including G-T wobbles
     
     
     for (nt=3; nt<=k-2; nt++) {  //% where is bulge, keep two bases on each side
      if ((S[i] != 0)&&(i < N)) {
	misWC = mWC & (-(1<<(2*nt)));
	//	printf(" %d", misWC);
	misWC /= 4;
	misWC += (mWC & ((1<<(2*nt-2))-1));
	misWC += (Swc[i]-1)*M/4;
	UpdateGTScores(C, p, P, misWC, k, kmin, imax, Ss, Cs, 0, 0, 1);
      }
     }

     // one-bulges, template-side, including G-T wobbles
     /*
     for (nt=3; nt<=k-2; nt++) {  //% where is bulge, keep two bases on each side
      for (bp = 0; bp < 4; bp++) {
	misWC = mWC & ((1<<(2*k-2))-(1<<(2*nt-2)));
       	misWC *= 4;
	misWC += (mWC & ((1<<(2*nt-2))-1));
	misWC += (bp<<(2*nt-2));
	UpdateGTScores(C, p, P, misWC, k, kmin, imax, Ss, Cs, 0, 0, 1);
      }
      }*/

     // two-bulges, sense-side, including G-T wobbles
         
     for (nt=4; nt<=k-2; nt++) {  //% where is bulge, keep two bases on each side
      if ((S[i] != 0)&&(S[i+1] != 0)&&(i < (N-1))) {
	misWC = mWC & (-(1<<(2*nt)));
	misWC /= 16;
	misWC += (mWC & ((1<<(2*nt-4))-1));
	misWC += (Swc[i]-1)*M/16;
	misWC += (Swc[i+1]-1)*M/4;
	UpdateGTScores(C, p, P, misWC, k, kmin, imax, Ss, Cs, 0, 0, 1);
      } 
     }

     // two-bulges, template-side, including G-T wobbles
     /*
     for (nt=3; nt<=k-3; nt++) {  //% where is bulge, keep two bases on each side
      for (bp = 0; bp < 4; bp++) {
       for (bp2 = 0; bp2 < 4; bp2++) {
	misWC = mWC & ((1<<(2*k-2))-(1<<(2*nt-2)));
       	misWC *= 16;
	misWC += (mWC & ((1<<(2*nt-2))-1));
	misWC += (bp<<(2*nt-2));
	misWC += (bp2<<(2*nt));
	UpdateGTScores(C, p, P, misWC, k, kmin, imax, Ss, Cs, 0, 0, 1);
       }
      }
     }
     */


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

      for (j=0; j<n; j++) {
        if (eq[j] == eq[i]) { S[j] = S[i]; marked[j]=1; }
      }

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

double verboten(char *Sorig)
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

 M = 1<<12;
 m = 0;                //% the ID of current subsequence
 for (i=0; i<N; i++) {
  if (S[i] == 0) {     //% reset strand & complex indicators
    m=0; K=5;          //% wait till full 6-mer considered
  } else {
    m = (  (m<<2) + S[i]-1 ) & (M-1);
    if (K) K=K-1; else V=V+bad6mer[m]; 
  }				       
 }

 free(S);
 return V; 
}

// --------------------------------------------------------------------------

// urn from Hofacker, Fontana 
unsigned short xsubi[3];  // for random number generator
double urn(void)    
                /* uniform random number generator; urn() is in [0,1] */
                /* uses a linear congruential library routine */ 
                /* 48 bit arithmetic */
{
    extern double erand48(unsigned short[3]);
    return erand48(xsubi);
}

int int_urn(int from, int to)
{
    return ( ( (int) (urn()*(to-from+1)) ) + from );
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

void mutate(char *S, char *St, int *wc, int *eq)
{
  int i; char oldc;
  do {
   i = int_urn(0,strlen(St)-1);
   oldc = S[i];
   S[i] = randbasec(St[i]);
   constrain(S,wc,eq);
  } while (oldc == S[i]);
  // if no base change allowed by template can alter constrained S, never halts
}
// --------------------------------------------------------------------------

double score_tmp(char *S, int *wc, int *eq)
{
  
  int factors[18] = {64, 512, 4096, 32768, 262144, 2097152,        // wuz powers of 5, now powers of 6, DZ 7/2/01
                      8,  64, 512,  4096,  32768, 262144,
                      1,   8,  64,  512,   4096,  32768 };
  int *C; double score=0.0; int i; double V;

  C = spurious(0,S,3,8,wc,eq);
  for (i=0; i<18; i++) score+=(factors[i]*C[i]*1.0)/3125.0;

  C = spurious(1,S,5,10,NULL,NULL);
  for (i=0; i<18; i++) score+=(factors[i]*C[i]*1.0)/3125.0;

  V = verboten(S);

  return score+V;
}
// --------------------------------------------------------------------------

void main(int argc, char *argv[])
{
 int debug=0;
 FILE *f; double r; int i,k,N=0,Nt=0; 
 int *testwc, *testeq; char *testS, *testSt; char c;
 int prob;
 int *C;
 time_t t;
 char *tempminS;
 int tempminscore;
 

 if (argc != 6) {
    fprintf(stderr,"%s rS St wc eq p\n",argv[0]);
    exit(-1);
 }

   time(&t);
   xsubi[0] = (unsigned short) t;
   xsubi[1] = (unsigned short) (t >> 16);
   xsubi[2] = 5246;

   prob = atoi(argv[5]);

 make_bad6mer();

 testS = malloc(10000); // big enough?
 if ( (f = fopen(argv[1],"r")) == NULL ) {
   fprintf(stderr,"Can't open rS file <%s>\n",argv[1]); exit(-1);
 } else {
   while ((c = fgetc(f)) != EOF) testS[N++]=c; testS[N]=0;
   if (debug) printf("S = <%s>\n",testS);
   N = strlen(testS);
   fclose(f);
 }

 testSt = malloc(10000);
 tempminS = malloc(10000);

 if ( (f = fopen(argv[2],"r")) == NULL ) {
   fprintf(stderr,"Can't open St file <%s>\n",argv[2]); exit(-1);
 } else {
   while ((c = fgetc(f)) != EOF) testSt[Nt++]=c; testSt[Nt]=0;
   if (debug) printf("St = <%s>\n",testSt);
   Nt = strlen(testSt);
   fclose(f);
 }

 if (N != Nt) 
   { fprintf(stderr,"S and St are different lengths!!!\n"); exit(-1); }

 testwc = calloc(N, sizeof(int));
 if ( (f = fopen(argv[3],"r")) == NULL ) {
   fprintf(stderr,"Can't open wc file <%s>\n",argv[3]); exit(-1);
 } else {
   for (i=0; i<N; i++) { fscanf(f," %lf",&r); testwc[i]=r; }
   if (debug) {
    printf("wc = <"); 
    for (i=0; i<N; i++) { printf("%d ", testwc[i]); }
    printf("> \n"); 
   }
   fclose(f);
 }

 testeq = calloc(N, sizeof(int));
 if ( (f = fopen(argv[4],"r")) == NULL ) {
   fprintf(stderr,"Can't open eq file <%s>\n",argv[4]); exit(-1);
 } else {
   for (i=0; i<N; i++) { fscanf(f,"%lf",&r); testeq[i]=r; }
   if (debug) {
    printf("eq = <"); 
    for (i=0; i<N; i++) { printf("%d ", testeq[i]); }
    printf("> \n"); 
   }
   fclose(f);
 }

 constrain(testS,testwc,testeq);
 printf("\nconstrained S = <\n%s\n>  N=%d \n",testS,N);

 C = spurious(0, testS, 3,8, testwc, testeq);
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

 { char *oldS; double old_score, new_score; int steps=0; double tempminscore = -1; int minstep=0;
   oldS = malloc(N); 
   old_score = score_tmp(testS,testwc,testeq);
   printf("%d steps: score = %11.5f \n",steps,old_score);
   for (i=0; i<N; i++) tempminS[i] = testS[i];
   tempminscore = old_score;
   minstep = steps;
  while (1) { 
   steps++;  if (steps%300 == 0) printf("\n\nStep%d  Current Best: %s\nScore:%f\nAt Step:%d\n\n",steps, tempminS, tempminscore, minstep);
   for (i=0; i<N; i++) oldS[i]=testS[i];
   mutate(testS,testSt,testwc,testeq);
   new_score = score_tmp(testS,testwc,testeq);
   if ((new_score <= old_score) || ((prob != 0) && (int_urn(0, prob) == 0))) {
     old_score = new_score;
     printf("%d steps: score = %11.5f \n",steps,new_score);
   } else {
     for (i=0; i<N; i++) testS[i]=oldS[i];
   }
   if (new_score < tempminscore) {
     minstep = steps;
     tempminscore = new_score;
     for (i=0; i<N; i++) tempminS[i]=testS[i];
   }
   if (steps - minstep > 100) { // eek... bad mutation... reload
     for (i=0; i<N; i++) testS[i]=tempminS[i];
     old_score = tempminscore;
     new_score = tempminscore;
   }
   if (steps - minstep > 2000) { // gee... stuck in a deep pit... stop program
     printf("\n\nFinal Answer: %s\nScore:%f\n\n", tempminS, tempminscore);
     exit(0);
   }
  }  
 }
}
