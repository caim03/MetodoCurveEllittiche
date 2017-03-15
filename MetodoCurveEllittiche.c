#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "modulus.h"

#define B1 1000 // B-smoothness bound

/* This struct define the structure of a point */
struct point
{
	long long x; // x-coordinate
	long long y; // y-coordinate
};

/* this function initialize the values needed for the algorithm */
void initialize(long long *n, struct point **P, struct point **Q,long long *A, long long *B, long *temp1, int *temp2){

	// n1 and n2 are prime numbers that define n to factorize
	long long n1, n2;						

	printf ("Insert a prime number n1: ");
	if (!scanf("%lld", &n1)){
		fprintf(stderr, "Error in scan\n");
		exit(EXIT_FAILURE);
	}

	printf ("Insert a prime number n2: ");
	if (!scanf ("%lld", &n2)){
		fprintf (stderr, "Error in scan\n");
		exit (EXIT_FAILURE);
	}

	// n is the number that should be factorized
	*n = (n1) * (n2);
	
	*P = malloc(sizeof(struct point));
	if (*P == NULL){
		fprintf (stderr, "Error in point initialize\n");
		exit (EXIT_FAILURE);
	}

	*Q = malloc(sizeof(struct point));
	if (*Q == NULL){
		fprintf (stderr, "Error in point initialize\n");
		exit (EXIT_FAILURE);
	}

	// Set random seed
	srand(time(NULL)); 

	(*P) -> x = rand() % *n;
	(*P) -> y = rand() % *n;

	*A = rand() % *n;
	*B = pow((*P) -> y, 2) - pow((*P) -> x, 3) - (*A)*((*P)->x);
	*B = (*B) % (*n);	

	int i;
	int j = 0;

	for (i = 2; i <= B1; i++){

		if (is_prime(i) == TRUE){
			// printf("%ld\n", factor);
			temp1[j] = i;
			// printf("temp1[j]=%ld\n", temp1[j]);
			j += 1;
			// printf("j = %d\n", j);
		}
	}

	for (i = 0; i < j; i++) {
		int m = 1;
		while (pow(temp1[i], m) <= B1)
			m++;
		
		temp2[i] = m - 1;
		// printf("temp2[%d] = %d\n", i, temp2[i]);
	}

	printf ("\n\nA = %lld, B = %lld, x0 = %lld, y0 = %lld\n\n", *A, *B, (*P) -> x, (*P) -> y);
}

short first_phase (int r, struct point *P, long long A, long long n, long long *p, 
	struct point *Q, long *pi, int *expo) {

	// indicate if a number has inverse or not
	long long a;
	// result of the first_phase	
	short res;		
	// value of mcd function between two numbers
	long long gcd;	
	// angular coefficient						
	long long m;	
	// index used in for cycle								
	int i;									

	Q -> x = P -> x;
	Q -> y = P -> y;
	
	for (i = 0; i < r; i++) {				
		long long k, power;
		power = pow(pi[i], expo[i]);
		for (k = 1; k <= power; k++){

			/* Case in which the points are egual */
			if (((Q -> x) == (P -> x)) && ((Q -> y) == (P -> y))) {	
				// find the inverse if exists	
				a = inverse(2 * (Q -> y), n); 

				/*  If the inverse not exists */
				if (a == -1) {							
					a = 2*(Q -> y);
					// modulus reduction
					a = a % n;
					gcd = mcd(a, n);
				
					if ((1 <= gcd) && (gcd <= n)) {
						// we have found a non-trivial factor of n
						*p = gcd;
						res = 1;
						return res;
					}
				}				

				// calculates the angular coefficient m = (3x^2 + A)(2y)^(-1)
				m = ((3 * pow(P -> x, 2) + A) * (a));
				// modulus reduction 
				m = m % n;		
			}

			else {									
				a = inverse((Q -> x) - (P -> x), n);

				/*  If the inverse not exists */
				if (a == -1) {						
					a = (Q -> x) - (P -> x);
					// modulus reduction
					a = a % n;
					gcd = mcd(a, n);
				
					if ((1 <= gcd) && (gcd <= n)) {
						// we have found a non-trivial factor of n
						*p = gcd;					
						res = 1;
						return res;
					}
					
				}

				// calculates the angular coefficient m = (y2-y1)*(x2-x1)
				m = mulmod((Q -> y) - (P -> y), a, n);
				m = m % n;
			}
				
			Q -> x = mulmod(m, m, n) - (Q -> x) - (P -> x);
			Q -> x = (Q -> x) % n;
			if ((Q -> x) < 0)
				Q -> x = n + (Q -> x);
			
			Q -> y = -(mulmod(m, (Q -> x) - (P -> x) + (P -> y), n));
			Q -> y = (Q -> y) % n;
			if ((Q -> y) < 0)
				Q -> y = n + (Q -> y);
		}	
	}

	res = 0;
	printf ("Non trivial factors are not found\n");
	return res;

}

int main(void){
	
	// Integer to be factored
	long long n, n1, n2;
	// Random point P in Zn x Zn 			
	struct point *P; 	
	// Integer of elliptic curve			
	long long A; 		
	// Integer of elliptic curve			
	long long B; 		
	// Discriminant of elliptic curve			
	long long delta; 				

	// index in for cycle
	int i;							
	long temp1[B1 * sizeof(long)];
	int temp2[B1 * sizeof(int)];

	// 0 or 1 to continue with second phase of to return p that is factor of n
	short res;			
	// is a factor of n, if the algorithm has success
	long long p;	
	// is a finite point that is used in second phase
	struct point *Q;	
	
	/* Initialize the main parameters of algorithm */
	initialize (&n, &P, &Q, &A, &B, temp1, temp2);
	int r = 0;
	while (temp1[r] != 0)
		r++;

	long pi[r];
	int expo[r];

	for (i = 0; i < r; i++) {
		pi[i] = temp1[i];
		expo[i] = temp2[i];
	}

	delta = 4 * pow(A, 3) + 27 * pow(B, 2);

	/* Check if the curve is an elliptic curve on Zn */
	if ((delta % n) == 0){
		printf ("This is not an elliptic curve on Z%lld\n", n);
		return 0;
	}

	printf("This is an elliptic curve on Z%lld\n", n);

	res = first_phase (r, P, A, n, &p, Q, pi, expo);

	if (res == 1){

		/* THE FIRST PHASE FINISHED CORRECTLY - RETURN p (factor of n) */
		printf("\t\t\t%lld is a non-trivial factor of %lld\n", p, n);
	}

	else {

		/* THE FIRST PHASE IS FAILED - START WITH SECOND PHASE */
		printf("A = %lld, B = %lld.\n", A, B);

		printf("Non avendo trovato fattori non banali relativi alla curva \nY^2 = X^3 + (%lld)*X + (%lld) e bound B1 %d, \n", A, B, B1);
		printf("la prima fase è fallita, procediamo con la seconda fase partendo\n dal punto al finito Q = (%lld, %lld).\n", (Q -> x), (Q -> y));
	}

	return 0;
}











