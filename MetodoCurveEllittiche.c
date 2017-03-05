#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define B1 1000 // B-smoothness bound

/* This struct define the structure of a point */
struct point
{
	long x; // x-coordinate
	long y; // y-coordinate
};

int is_prime (int n) {	

	/* this function take in input a number and check if it is prime or not */

	int i;

	if (n == 1 || n == 2)
			return 1;
		

	for (i = 2; i < n; i++){

		if (n % i == 0)
			return 0;
		
	}

	return n;
}

void initialize(long *n, struct point **P, long *A, long *B, long *temp1, int *temp2){
	/* this function initialize the values needed for the algorithm */

	long n1, n2;						// n1 and n2 are prime numbers that define n to factorize
	long factor;

	printf ("Inserire numero primo n1: ");
	if (!scanf("%ld", &n1)){
		fprintf(stderr, "Error in scan\n");
		exit(EXIT_FAILURE);
	}

	printf ("Inserire numero primo n2: ");
	if (!scanf ("%ld", &n2)){
		fprintf (stderr, "Error in scan\n");
		exit (EXIT_FAILURE);
	}

	*n = (n1) * (n2);			// n is the number that should be factorized with algorithm
	
	*P = malloc (sizeof(struct point));
	if (*P == NULL){
		fprintf (stderr, "Error in point initialize\n");
		exit (EXIT_FAILURE);
	}

	(*P) -> x = rand() % *n;
	(*P) -> y = rand() % *n;

	*A = rand() % *n;
	*B = pow((*P) -> y, 2) - pow((*P) -> x, 3) - (*A)*((*P)->x);
	*B = (*B) % (*n);	

	int i;
	int j = 0;

	for (i = 2; i <= B1; i++){

		factor = is_prime(i);
		if ((factor != 0) && (factor != 1)){
			printf("%ld\n", factor);
			temp1[j] = factor;
			printf("temp1[j]=%ld\n", temp1[j]);
			j += 1;
			printf("j = %d\n", j);
		}

	}

	for (i = 0; i < j; i++) {
		int m = 1;
		while (pow(temp1[i], m) <= B1)
			m++;
		
		temp2[i] = m - 1;
		printf("temp2[%d] = %d\n", i, temp2[i]);
	}

	printf ("\n\nA = %ld, B = %ld, x0 = %ld, y0 = %ld\n\n", *A, *B, (*P) -> x, (*P) -> y);
}


long mcd(long a, long b) {
	/* this function calculate gcd between number a and number b */

	if (a < 0)
		a = b + a;

	return (b != 0)?mcd(b, a % b) : a;
}


long inverse (long a, long b) {
	/* this function calculate inverse of a modulus b */

	long b0 = b, t, q;
	long x0 = 0, x1 = 1;

	a = a % b;
	if (a < 0)
		a = b + a;

	if (mcd(a, b) != 1)
		return -1;

	while (a > 1) {
		q = a / b;
		t = b, b = a % b, a = t;
		t
		 = x0, x0 = x1 - q * x0, x1 = t;
	}
	if (x1 < 0) 
		x1 += b0;

	return x1;
}


short first_phase (long k, struct point *P, long A, long n, long *p, struct point *Q) {
	
	/* da quello che ho capito, facciamo (k * P) e ogni volta controlliamo se gcd (P -> y, n)!= 1;
	   se è così allora continuiamo fino a fare (k * P) e non trovando nessuna P -> y invertibile
	   avremo Q = (k * P) che useremo nella seconda fase. Se facendo queste moltiplicazioni troviamo
	   una 2*(Q -> y) oppure ((Q -> x) - (P -> x)) che non ha inverso allora il gcd sarà il fattore 
	   non banale trovato, quindi l'algoritmo ha successo e ritorniamo il gcd trovato. */

	long a;									// indicate if a number has inverse or not
	short res;								// result of the first_phase
	long gcd;								// value of mcd function between two numbers
	long m;									// angular coefficient
	int i;									// index used in for cycle

	printf("P = (%ld, %ld)\n", P -> x, P -> y);


	Q -> x = P -> x;
	Q -> y = P -> y;

	printf("Q = (%ld, %ld)\n\n", Q -> x, Q -> y);

	
	for (i = 1; i <= k; i++) {				

		if (i == 1) {		

			a = inverse (2 * (Q -> y), n); // calcolo l'inverso se esiste

			/*  Se non trovo l'inverso */
			if (a == -1) {	
				printf ("%d. %ld non ha inverso mod %ld nell'if\n", i, 2 * (Q -> y), n);
					
				a = 2*(Q -> y);
				a = a % n;
				gcd = mcd(a, n);

				printf ("%d gcd (%ld, %ld) = %ld nell'if\n", i, a, n, gcd);
			
				if ((1 <= gcd) && (gcd <= n)) {
					*p = gcd; // fattore non banale di n
					res = 1;
					return res;	// algoritmo ha avuto successo, maria io esco.
				}

			printf ("%d. Inverso di %ld mod %ld è a = %ld nell'if\n", i, 2 * (Q -> y), n, a);
				
			}

			m = ((3 * pow(P -> x, 2) + A) * (a)); // calcolo m = (3x^2 + A)(2y)^(-1)
			m = m % n;		
		}

		else {									// caso in cui cambia m
			/* qui devo farel'inverso di (x2-x1) quindi mi muovo diversamente*/
			a = inverse ((Q -> x) - (P -> x), n);
			printf ("%d. Inverso di %ld mod %ld è a = %ld nell'else\n", i, (Q -> x) - (P -> x), n, a);

			/* se non trovo l'inverso */
			if (a == -1) {						
				a = (Q -> x) - (P -> x);
				a = a % n;
				gcd = mcd(a, n);
				printf ("%d gcd (%ld, %ld) = %ld nell'else\n", i, a, n, gcd);
			
				if ((1 <= gcd) && (gcd <= n)) {
					*p = gcd;					// fattore non banale di n
					printf ("%d. Abbiamo trovato un fattore non banale con Q = (%ld, %ld)\n\n", i, Q -> x, Q -> y);
					res = 1;
					return res;					// algoritmo ha avuto successo, maria io esco.
				}
				
			}

			m = ((Q -> y) - (P -> y)) * a; 		//calcolo m = (y2-y1)*(x2-x1)
			m = m % n;
		}
		
		Q -> x = pow(m, 2) - (Q -> x) - (P -> x);
		Q -> x = (Q -> x) % n;
		if ((Q -> x) < 0)
			Q -> x = n + (Q -> x);
		

		Q -> y = -(m * ((Q -> x) - (P -> x)) + (P -> y));
		Q -> y = (Q -> y) % n;
		if ((Q -> y) < 0)
			Q -> y = n + (Q -> y);
		

		printf ("Q = (%ld, %ld)\n", Q -> x, Q -> y);
	}

	res = 0;
	printf ("Non ho trovato fattori non banali.\n");
	return res;

}


long kvalue (long *pi, int *expo, long r) {

	/* this function calculate the value of k, that indicates the number of time that should sum
	   a point P in first phase of the algorithm. */

	//int dim;									// dimension of the array pi/expo
	int i;										// index in for cycle
	long k = 1;									


	//dim = ( sizeof(pi) / sizeof(long)) + 1;
	printf ("Il numero di fattori primi è %ld.\n\n", r);

	for (i = 0; i < r; i++) {
		k = k * pow(pi[i], expo[i]);
		printf("%d k = %ld\n", i, k);
	}	
	
	
	printf ("Si effettueranno al massimo k = %ld somme di P alla ricerca di fattori non banali.\n\n", k);
	return k;
}


int main(void){
	
	long n, n1, n2; 						// Integer to be factored
	struct point *P; 						// Random point P in Zn x Zn
	long A; 								// Integer of elliptic curve
	long B; 								// Integer of elliptic curve
	long k; 								// Large B-smooth number
	long delta; 							// Discriminant of elliptic curve
	//int *expo;								// array of exponents for prime numbers of pi
	//long *pi;								// array with prime numbers <= BSB
	int i;									// index in for cycle
	long temp1[B1*sizeof(long)];
	int temp2[B1*sizeof(int)];

	//printf("B1*sizeof(long) = %ld\n", B1*sizeof(long));
	




	/* Initialize the main parameters of algorithm */
	initialize (&n, &P, &A, &B, temp1, temp2);
	long r = 0;
	while (temp1[r] != 0)
		r++;
	printf("r = %ld\n", r);
	//pi= malloc (r*sizeof(long));		
	//expo = malloc (r*sizeof(int));
	long pi[r];
	int expo[r];

	//printf("dimensione pi = %ld\n", sizeof(&pi));

	for (i = 0; i < r; i++) {
		pi[i] = temp1[i];
		expo[i] = temp2[i];

		printf ("Fattori primi pi[%d] = %ld,	con esponente expo[%d] = %d\n", i, pi[i], i, expo[i]);
	}

	k = kvalue (pi, expo, r);
	delta = 4 * pow(A, 3) + 27 * pow(B, 2);

	/* Check if the curve is an elliptic curve on Zn */
	if ((delta % n) == 0){
		printf ("This is not an elliptic curve on Z%ld\n", n);
		return 0;
	}

	printf ("This is an elliptic curve on Z%ld\n", n);

	/* Start first phase */
	short res;			// 0 or 1 to continue with second phase of to return p that is factor of n
	long p;				// is a factor of n, if the algorithm has success
	struct point *Q;	// is a finite point that is used in second phase

	Q = malloc(sizeof(struct point));
	if (Q == NULL){
		fprintf (stderr, "Error in point initialize\n");
		exit (EXIT_FAILURE);
	}

	res = first_phase (k, P, A, n, &p, Q);

	if (res == 1){

		/* THE FIRST PHASE FINISHED CORRECTLY - RETURN p (factor of n) */
		printf ("\t\t\t%ld è un fattore non banale di %ld\n", p, n);

	}

	else {

		/* THE FIRST PHASE IS FAILED - START WITH SECOND PHASE */
		printf ("A = %ld, B = %ld.\n", A, B);

		printf ("Non avendo trovato fattori non banali relativi alla curva \nY^2 = X^3 + (%ld)*X + (%ld) e bound B1 %d, \n", A, B, B1);
		printf ("la prima fase è fallita, procediamo con la seconda fase partendo\n dal punto al finito Q = (%ld, %ld).\n", (Q -> x), (Q -> y));

	}

	return 0;
}











