#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define B1 1000 // B-smoothness bound

/* This struct define the structure of a point */
struct point
{
	long long x; // x-coordinate
	long long y; // y-coordinate
};

/* This function performs a modulus moltiplication to prevent overflow problems */
long long mulmod(long long a, long long b, long long m){
	long long r = 0;

	while (b > 0) {
        if (b & 1)  r = ((m-r) > a) ? r+a : r+a-m;    /* r = (r + a) % m */
        b >>= 1;
        if (b)      a = ((m-a) > a) ? a+a : a+a-m;    /* a = (a + a) % m */
    }
    return r;
}

int is_prime (int n) {	

	/* this function take in input a number and check if it is prime or not */

	int i;

	if (n == 1 || n == 2)
			return n;
		

	for (i = 2; i < n/2; i++){

		if (n % i == 0)
			return 0;
		
	}

	return n;
}

void initialize(long long *n, struct point **P, long long *A, long long *B, long *temp1, int *temp2){
	/* this function initialize the values needed for the algorithm */

	long long n1, n2;						// n1 and n2 are prime numbers that define n to factorize
	long factor;

	printf ("Inserire numero primo n1: ");
	if (!scanf("%lld", &n1)){
		fprintf(stderr, "Error in scan\n");
		exit(EXIT_FAILURE);
	}

	printf ("Inserire numero primo n2: ");
	if (!scanf ("%lld", &n2)){
		fprintf (stderr, "Error in scan\n");
		exit (EXIT_FAILURE);
	}

	*n = (n1) * (n2);			// n is the number that should be factorized with algorithm
	
	*P = malloc (sizeof(struct point));
	if (*P == NULL){
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

	printf ("\n\nA = %lld, B = %lld, x0 = %lld, y0 = %lld\n\n", *A, *B, (*P) -> x, (*P) -> y);
}


long mcd(long long a, long long b) {
	/* this function calculate gcd between number a and number b */
	if (a == 0){
		return b;
	}

	if (b == 0){
		return a;
	}

	if (a < 0)
		a = b + a;

	return (b != 0)?mcd(b, a % b) : a;
}


long long inverse (long long a, long long b) {
	/* this function calculate inverse of a modulus b */

	long long b0 = b, t, q;
	long long x0 = 0, x1 = 1;

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


short first_phase (int r, struct point *P, long long A, long long n, long long *p, 
	struct point *Q, long *pi, int *expo) {
	
	/* da quello che ho capito, facciamo (k * P) e ogni volta controlliamo se gcd (P -> y, n)!= 1;
	   se è così allora continuiamo fino a fare (k * P) e non trovando nessuna P -> y invertibile
	   avremo Q = (k * P) che useremo nella seconda fase. Se facendo queste moltiplicazioni troviamo
	   una 2*(Q -> y) oppure ((Q -> x) - (P -> x)) che non ha inverso allora il gcd sarà il fattore 
	   non banale trovato, quindi l'algoritmo ha successo e ritorniamo il gcd trovato. */

	long long a;							// indicate if a number has inverse or not
	short res;								// result of the first_phase
	long long gcd;							// value of mcd function between two numbers
	long long m;									// angular coefficient
	int i;									// index used in for cycle

	printf("P = (%lld, %lld)\n", P -> x, P -> y);


	Q -> x = P -> x;
	Q -> y = P -> y;

	printf("Q = (%lld, %lld)\n\n", Q -> x, Q -> y);

	
	for (i = 0; i < r; i++) {				
		long long k, power;
		power = pow(pi[i], expo[i]);
		for (k = 1; k <= power; k++){
			if (((Q -> x) == (P -> x)) && ((Q -> y) == (P -> y))) {		
				a = inverse (2 * (Q -> y), n); // calcolo l'inverso se esiste

				/*  Se non trovo l'inverso */
				if (a == -1) {	
					printf ("%lld non ha inverso mod %lld nell'if\n", 2 * (Q -> y), n);
						
					a = 2*(Q -> y);
					a = a % n;
					gcd = mcd(a, n);

					printf ("gcd (%lld, %lld) = %lld nell'if\n",a, n, gcd);
				
					if ((1 <= gcd) && (gcd <= n)) {
						*p = gcd; // fattore non banale di n
						res = 1;
						return res;	// algoritmo ha avuto successo, maria io esco.
					}
				}

				printf ("Inverso di %lld mod %lld è a = %lld nell'if\n", 2 * (Q -> y), n, a);
				

				m = ((3 * pow(P -> x, 2) + A) * (a)); // calcolo m = (3x^2 + A)(2y)^(-1)
				m = m % n;		
			}

			else {									
				/* qui devo fare l'inverso di (x2-x1) quindi mi muovo diversamente*/
				a = inverse ((Q -> x) - (P -> x), n);
				printf ("Inverso di %lld mod %lld è a = %lld nell'else\n", (Q -> x) - (P -> x), n, a);

				/* se non trovo l'inverso */
				if (a == -1) {						
					a = (Q -> x) - (P -> x);
					a = a % n;
					gcd = mcd(a, n);
					printf ("gcd (%lld, %lld) = %lld nell'else\n", a, n, gcd);
				
					if ((1 <= gcd) && (gcd <= n)) {
						*p = gcd;					// fattore non banale di n
						printf ("Abbiamo trovato un fattore non banale con Q = (%lld, %lld)\n\n", Q -> x, Q -> y);
						res = 1;
						return res;					// algoritmo ha avuto successo, maria io esco.
					}
					
				}

				m = mulmod((Q -> y) - (P -> y), a, n); 		//calcolo m = (y2-y1)*(x2-x1)
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
			

			printf ("Q = (%lld, %lld)\n", Q -> x, Q -> y);
		}	
	}

	res = 0;
	printf ("Non ho trovato fattori non banali.\n");
	return res;

}


int main(void){
	
	long long n, n1, n2; 			// Integer to be factored
	struct point *P; 				// Random point P in Zn x Zn
	long long A; 					// Integer of elliptic curve
	long long B; 					// Integer of elliptic curve
	long long delta; 				// Discriminant of elliptic curve

	int i;							// index in for cycle
	long temp1[B1*sizeof(long)];
	int temp2[B1*sizeof(int)];
	
	/* Initialize the main parameters of algorithm */
	initialize (&n, &P, &A, &B, temp1, temp2);
	int r = 0;
	while (temp1[r] != 0)
		r++;
	printf("r = %d\n", r);

	long pi[r];
	int expo[r];

	for (i = 0; i < r; i++) {
		pi[i] = temp1[i];
		expo[i] = temp2[i];

		printf ("Fattori primi pi[%d] = %ld,	con esponente expo[%d] = %d\n", i, pi[i], i, expo[i]);
	}

	delta = 4 * pow(A, 3) + 27 * pow(B, 2);

	/* Check if the curve is an elliptic curve on Zn */
	if ((delta % n) == 0){
		printf ("This is not an elliptic curve on Z%lld\n", n);
		return 0;
	}

	printf ("This is an elliptic curve on Z%lld\n", n);

	/* Start first phase */
	short res;			// 0 or 1 to continue with second phase of to return p that is factor of n
	long long p;	// is a factor of n, if the algorithm has success
	struct point *Q;	// is a finite point that is used in second phase

	Q = malloc(sizeof(struct point));
	if (Q == NULL){
		fprintf (stderr, "Error in point initialize\n");
		exit (EXIT_FAILURE);
	}

	res = first_phase (r, P, A, n, &p, Q, pi, expo);

	if (res == 1){

		/* THE FIRST PHASE FINISHED CORRECTLY - RETURN p (factor of n) */
		printf ("\t\t\t%lld è un fattore non banale di %lld\n", p, n);

	}

	else {

		/* THE FIRST PHASE IS FAILED - START WITH SECOND PHASE */
		printf ("A = %lld, B = %lld.\n", A, B);

		printf ("Non avendo trovato fattori non banali relativi alla curva \nY^2 = X^3 + (%lld)*X + (%lld) e bound B1 %d, \n", A, B, B1);
		printf ("la prima fase è fallita, procediamo con la seconda fase partendo\n dal punto al finito Q = (%lld, %lld).\n", (Q -> x), (Q -> y));

	}

	return 0;
}











