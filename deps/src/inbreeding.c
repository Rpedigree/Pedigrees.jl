#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>

/**
 * Create the inbreeding coefficients according to the algorithm given
 * in "Comparison of four direct algorithms for computing inbreeding
 * coefficients" by Mehdi Sargolzaei and Hiroaki Iwaisaki, Animal
 * Science Journal (2005) 76, 401--406.  This function is a modified
 * version of the code published in an appendix to that paper.
 *
 * @param sire 1-based indices of animal's sire (0 => unknown)
 * @param dam  1-based indices of animal's dam (0 => unknown)
 * @param lap  longest ancestral path.  lap[0] == -1, for unknown parent
 * @param n    length of sire and dam.  lap and F are of length n+1
 * @param F    double vector of length n that contains the inbreeding coefficients on return
 *
 */

void printivec(const ptrdiff_t* vv, const char* nm, ptrdiff_t l) {
    ptrdiff_t i;
    printf("%s = [%ld", nm, vv[0]);
    for(i = 1; i < l; ++i) printf(", %ld", vv[i]);
    printf("]\n");
}

void inbreeding(const ptrdiff_t* sire, const ptrdiff_t* dam, const ptrdiff_t *lap,
		ptrdiff_t n, double *F) {
    ptrdiff_t i, j, t, maxlap, S, D;
    ptrdiff_t *SI, *MI;		/* start and minor */
    double *L = calloc(n + 1, sizeof(double)), *B = calloc(n + 1, sizeof(double));
    ptrdiff_t *Anc = calloc(n + 1, sizeof(ptrdiff_t));	/* ancestor */
	
    maxlap = lap[0];	     /* determine the maximum ancestral path length */
    for (i = 1; i <= n; ++i) {
	if (lap[i] > maxlap) maxlap = lap[i];
    }
    SI = calloc(maxlap + 1, sizeof(ptrdiff_t)); /* initialized to zeros */
    MI = calloc(maxlap + 1, sizeof(ptrdiff_t));

    F[0] = -1.0;		     /* set F  for unknown parents */
    for(i = 1; i <= n; i++) {	     /* evaluate F */
	S = sire[i-1]; D = dam[i-1]; /* parents of animal i */
	printf("(i,S,D) => (%ld,%ld,%ld)\n",i,S,D);
	B[i] = 0.5 - 0.25 * (F[S] + F[D]); 
				/* adjust start and minor */
	for (j = 0; j < lap[i]; j++) {++SI[j]; ++MI[j];} 
	printivec(SI,"SI",maxlap);
	printivec(MI,"MI",maxlap);
	if (S == 0 || D == 0) { /* both parents unknown */
	    printf("unknown parents exit\n");
	    F[i] = L[i] = 0.0; 
	    continue;
	}
	if(S == sire[i-2] && D == dam[i-2]) { /* full-sib with last animal */
	    printf("full sib exit\n");
	    F[i] = F[i-1];
	    L[i] = L[i-1];
	    continue;
	}
    	F[i] = -1.0; L[i] = 1.0; 
	t = lap[i]; /* largest lap group number in the animal's pedigree */
	Anc[MI[t]++] = i; /* initialize Anc and increment MI[t] */
	while(t > -1) { /* from the largest lap group to zero */
	    j = Anc[--MI[t]]; /* next ancestor */
	    printf("(i,j,t) => (%ld,%ld,%ld)\n",i,j,t);
	    printivec(SI,"SI",maxlap);
	    printivec(MI,"MI",maxlap);
	    printivec(Anc,"Anc",n+1);
	    S = sire[j-1]; D = dam[j-1]; /* parents of the ancestor */
	    if (S) {
		if (!L[S]) Anc[MI[lap[S]]++] = S; 
				/* add sire in its lap group in Anc
				 * array if it is not added yet and
				 * increment the minor index for the group */ 
		L[S] += 0.5 * L[j]; /* contribution to sire */
	    }
	    if (D) {
		if (!L[D]) Anc[MI[lap[D]]++] = D;
		L[D] += 0.5 * L[j]; /* contribution to dam */
	    }
	    F[i] += L[j] * L[j] * B[j];
	    L[j] = 0; /*clear L[j] for the evaluation of the next animal */
	    if (MI[t] == SI[t]) --t; /* move to the next lap group when
				      * all ancestors in group t have been
				      * evaluated */
	} 
    }
    free(L); free(B); free(Anc); free(SI); free(MI);
}
