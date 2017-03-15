#include <stdio.h>
#include <stdlib.h>
#include "modulus.h"

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

/* this function take in input a number and check if it is prime or not */
short is_prime (int n) {	
	int i;

	if (n == 1 || n == 2)
			return TRUE;

	for (i = 2; i < n/2; i++){
		// if i divides n
		if (n % i == 0)
			return FALSE;
		
	}

	return TRUE;
}

/* this function calculate gcd between number a and number b */
long long mcd(long long a, long long b) {
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

/* this function calculate inverse of a modulus b */
long long inverse (long long a, long long b) {

	long long b0 = b, t, q;
	long long x0 = 0, x1 = 1;

	// modulus reduction
	a = a % b;
	if (a < 0)
		a = b + a;

	// if inverse not exists
	if (mcd(a, b) != 1)
		return -1;

	while (a > 1) {
		q = a / b;
		t = b, b = a % b, a = t;
		t = x0, x0 = x1 - q * x0, x1 = t;
	}
	if (x1 < 0) 
		x1 += b0;

	return x1;
}