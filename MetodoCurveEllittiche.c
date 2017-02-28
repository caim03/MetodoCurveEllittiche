#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define BSB 3 // B-smoothness bound

/* This struct define the structure of a point */
struct point
{
	long x; // x-coordinate
	long y; // y-coordinate
};

void initialize(long *n, struct point *P, long *A, long *B){

	if(!scanf("%ld", n)){
		fprintf(stderr, "Error in scan\n");
		exit(EXIT_FAILURE);
	}
	
	P = malloc(sizeof(struct point));
	if (P == NULL){
		fprintf(stderr, "Error in point initialize\n");
		exit(EXIT_FAILURE);
	}

	P->x = rand() % *n;
	P->y = rand() % *n;

	*A = rand() % *n;
	*B = pow(P->y, 2) - pow(P->x, 3) - (*A)*(P->x);
	*B = *B % *n;	

	printf("A = %ld, B = %ld, x0 = %ld, y0 = %ld\n", *A, *B, P->x, P->y);
}

long factorial(void){
	long k = BSB;
	int i;

	for(i = 1; i < BSB; i++){
		k = k*(BSB-i);
	}

	return k;
}

int main(void){
	
	long n; // Integer to be factored
	struct point *P; // Random point P in Zn x Zn
	long A; // Integer of elliptic curve
	long B; // Integer of elliptic curve
	long k; // Large B-smooth number
	long delta; // Discriminant of elliptic curve

	/* Initialize the main parameters of algorithm */
	initialize(&n, P, &A, &B);
	k = factorial();
	delta = 4 * pow(A, 3) + 27 * pow(B, 2);

	/* Check if the curve is an elliptic curve on Zn */
	if ((delta % n) == 0){
		printf("This is not an elliptic curve on Z%ld\n", n);
		return 0;
	}

	printf("This is an elliptic curve on Z%ld\n", n);

	/* Start first phase */
	short res;
	long p;
	struct point *Q;

	Q = malloc(sizeof(struct point));
	if (Q == NULL){
		fprintf(stderr, "Error in point initialize\n");
		exit(EXIT_FAILURE);
	}

	res = first_phase(k, P, A, n, *p, Q);

	if (res == 1){

		/* THE FIRST PHASE FINISHED CORRECTLY - RETURN p (factor of n) */

	}

	else {

		/* THE FIRST PHASE IS FAILED - START WITH SECOND PHASE */

	}

	return 0;
}