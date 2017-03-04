#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define BSB 1000 // B-smoothness bound

/* This struct define the structure of a point */
struct point
{
	long x; // x-coordinate
	long y; // y-coordinate
};

int is_prime (int n) {		
	int i;

	if (n == 1 || n == 2){
			return 1;
		}

	for (i = 2; i < n; i++){

		if (n%i == 0){
			return 0;
		}
	}
	return 1;
}

void initialize(long *n, struct point **P, long *A, long *B, long *pi, int *expo){
	/* this function initialize the values needed for the algorithm*/

	long n1, n2;
	printf("Inserire numero primo n1: ");
	if(!scanf("%ld", &n1)){
		fprintf(stderr, "Error in scan\n");
		exit(EXIT_FAILURE);
	}

	printf("Inserire numero primo n2: ");
	if(!scanf("%ld", &n2)){
		fprintf(stderr, "Error in scan\n");
		exit(EXIT_FAILURE);
	}

	*n = (n1)*(n2);
	
	*P = malloc(sizeof(struct point));
	if (*P == NULL){
		fprintf(stderr, "Error in point initialize\n");
		exit(EXIT_FAILURE);
	}

	(*P)->x = rand() % *n;
	(*P)->y = rand() % *n;

	*A = rand() % *n;
	*B = pow((*P)->y, 2) - pow((*P)->x, 3) - (*A)*((*P)->x);
	*B = *B % *n;	

	int i;
	int j = 0;

	for (i=2; i <= BSB; i++){

		if (is_prime(i) == 1){
			pi[j] = i;
			j +=1;
		}

	}

	for (i = 0; i < j; i++) {
		int m = 1;
		while(pow(pi[i], m) <= BSB){
			m++;
		}
		expo[i] = m - 1;
	}

	printf("A = %ld, B = %ld, x0 = %ld, y0 = %ld\n", *A, *B, (*P)->x, (*P)->y);
}


long mcd(long a, long b) {

  long c;
  if (a < 0)
  	a = -a;

  if (b < 0)
  	b = -b;

  while (a != 0) {

  	c = a;
  	a = b % a;
  	b = c;

  }
  printf("b = %ld\n", b);
  return b;
}


long inverse (long a, long n) {
	long inv;							//qua salvo inverso di a
	long i;								//indice del for
	long cong;							//serve per trasformare numero negativo in positivo mod n

	if (a < 0)							//controllo se a è negativo per assegnare a cong il suo congruo mod n
		cong = n - a;
	else
		cong = a;


	for (i = 1; i < n; i++) {			//cerco il primo i che moltiplicato per a mi dia resto 1 mod n
		if (((cong*i) % n) == 1) {
			inv = i;					//ho trovato l'inverso cercato mod n
			printf("inverso di n %ld\n", inv);
			return inv;					
		}
	}

	return -1;

}



short first_phase (long k, struct point *P, long A, long n, long *p, struct point *Q) {
	/* da quello che ho capito, facciamo k*P e ogni volta controlliamo se gcd (P->y,n)!= 1;
	   se è così allora continuiamo fino a fare k*P e non trovando nessuna P->y invertibile
	   avremo Q = k*P che useremo nella seconda fase. Se facendo queste moltiplicazioni troviamo
	   una P->y che non ha inverso allora il gcd sarà il fattore non banale trovato, quindi
	   l'algoritmo ha successo e ritorniamo il gcd trovato. */
	long a;									// serve per calcolare ogni volta, se esiste, l'inverso
	short res;								// è la variabile che ritorniamo
	long gcd;
	long m;								// calcolo il fattore non banale e lo salvo in p
	int i;			

	printf("P -> x = %ld\n", P -> x);
	printf("P -> y = %ld\n", P -> y);


	Q -> x = P -> x;
	Q -> y = P -> y;

	printf("Q -> x = %ld\n", Q -> x);
	printf("Q -> y = %ld\n", Q -> y);
	
	for (i=1; i <= k; i++) {				

		if (i == 1) {		

			a = inverse (2*(Q -> y), n);			// la funzione che controlla se esiste inverso
			printf("inverso a = %ld nell'if\n", a);

			if (a == -1) {						// non ho trovato inverso modulo n
				a = 2*(Q -> y);
				a = a % n;
				gcd = mcd(a, n);

				printf("mcd relativo ad a = %ld nell'if\n", gcd);
			
				if ((1 <= gcd) && (gcd <= n)) {
					*p = gcd;					// fattore non banale di n
					res = 1;
					return res;					// algoritmo ha avuto successo, maria io esco.
				}

				
			}

			m = ((3*pow(Q -> x, 2) + A) * (a));
			m = m % n;

			Q -> x = pow(m, 2) - 2*(Q -> x);			// in Q salvo ogni volta quello che calcolo
			Q -> x = (Q -> x) % n;

			Q -> y = (-m)*(((Q -> x) - (P -> x)) + (P -> y));
			Q -> y = (Q -> y) % n;

			printf("Q -> x = %ld nell'if\n", Q -> x);
			printf("Q -> y = %ld nell'if\n\n", Q -> y);			
			
		}

		else {									// caso in cui cambia m
			/* qui devo farel'inverso di (x2-x1) quindi mi muovo diversamente*/
			a = inverse ((Q -> x) - (P -> x), n);
			printf("inverso a = %ld nell'else\n", a);


			if (a == -1) {						// non ho trovato inverso modulo n
				a = (Q -> x) - (P -> x);
				a = a % n;
				gcd = mcd(a, n);
				printf("gcd = %ld nell'else\n", gcd);
			
				if ((1 <= gcd) && (gcd <= n)) {
					*p = gcd;					// fattore non banale di n
					printf("Abbiamo trovato un fattore non banale con Q = (%ld, %ld)\n\n", Q -> x, Q -> y);
					res = 1;
					return res;					// algoritmo ha avuto successo, maria io esco.
				}
				
			}

			m = ((Q -> y) - (P -> y)) * a;		//calcolo m = (y2-y1)*(x2-x1)
			m = m % n;

			Q -> x = pow(m, 2) - (Q -> x) - (P -> x);
			Q -> x = (Q -> x ) % n;

			Q -> y = -(m * ((Q -> x) - (P -> x)) + (P -> y));
			Q -> y = (Q -> y ) % n;

			printf("Q -> x = %ld nell'else\n", Q -> x);
			printf("Q -> y = %ld nell'else\n\n", Q -> y);	

		}

		/*se arrivo qui significa che avrò un punto Q al finito che ho salvato come struttura e che userò
		  nelle seconda fase */
		

	}

	res = 0;
	return res;

}


long kvalue(long *pi, int *expo) {
	int dim;
	int i;
	long k = 1;

	dim = (sizeof(pi) / sizeof (long)) + 1;
	printf("dim pi & expo = %d\n", dim);

	for (i = 0; i < dim; i++){
		k = k * pow(pi[i], expo[i]);
	}
	
	printf("k = %ld\n", k);
	return k;
}


int main(void){
	
	long n, n1, n2; // Integer to be factored
	struct point *P; // Random point P in Zn x Zn
	long A; // Integer of elliptic curve
	long B; // Integer of elliptic curve
	long k; // Large B-smooth number
	long delta; // Discriminant of elliptic curve
	int *expo;
	long *pi;

	pi= malloc (BSB * sizeof(long));		
	expo = malloc (BSB * sizeof(int));


	/* Initialize the main parameters of algorithm */
	initialize(&n, &P, &A, &B, pi, expo);
	printf("P -> x = %ld\n", P -> x);

	int i;
	for (i=0; i <= (sizeof(pi)/sizeof(long)); i++){
		printf("pi[%d] = %ld,	expo[%d] = %d\n", i, pi[i], i, expo[i]);
	}

	k = kvalue(pi, expo);
	delta = 4 * pow(A, 3) + 27 * pow(B, 2);

	/* Check if the curve is an elliptic curve on Zn */
	if ((delta % n) == 0){
		printf("This is not an elliptic curve on Z%ld\n", n);
		return 0;
	}

	printf("This is an elliptic curve on Z%ld\n", n);

	/* Start first phase */
	short res;	// 0 or 1 to continue with second phase of to return p that is factor of n
	long p;		// is a factor of n, if the algorithm has success
	struct point *Q;	// is a finite point that is used in second phase

	Q = malloc(sizeof(struct point));
	if (Q == NULL){
		fprintf(stderr, "Error in point initialize\n");
		exit(EXIT_FAILURE);
	}

	res = first_phase(k, P, A, n, &p, Q);

	if (res == 1){

		/* THE FIRST PHASE FINISHED CORRECTLY - RETURN p (factor of n) */
		printf ("p = %ld è un fattore non banale di %ld\n", p, n);

	}

	else {

		/* THE FIRST PHASE IS FAILED - START WITH SECOND PHASE */

	}

	return 0;
}











