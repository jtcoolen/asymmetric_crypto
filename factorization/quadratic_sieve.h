#ifndef QUADRATIC_SIEVE_H
#define QUADRATIC_SIEVE_H

#include <gmp.h>
#include <stdint.h>

int64_t shanks_tonelli(int64_t a, int64_t p);

void gauss_eliminate(char *a, char *b, char *x, int n);

void quadratic_sieve(mpz_t n, unsigned long long P, unsigned long long A);

#endif
