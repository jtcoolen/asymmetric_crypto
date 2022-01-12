#include <gmp.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>

void f(mpz_t r, mpz_t a, mpz_t n) {
  mpz_mul(r, a, a);
  mpz_add_ui(r, r, 1);
  mpz_mod(r, r, n);
}

void factor_rho_pollard(mpz_t n, mpz_t factor)
{
  mpz_t xk, xj, diff, gcd;
  mpz_inits(xk, xj, diff, gcd, NULL);

  mpz_set_ui(xk, 1);
  mpz_set_ui(xj, 1);

  uint64_t k = 2;

  for(uint64_t j = 0;j < UINT64_MAX;j++) {
    f(xj, xj, n);

    if (j == k) {
      mpz_set(xk, xj);
      k <<= 1;
    }

    mpz_sub(diff, xk, xj);
    mpz_gcd(gcd, diff, n);
    if (mpz_cmp_ui(gcd, 1) != 0 && mpz_cmp(gcd, n) != 0) {
      mpz_set(factor, gcd);
      return;
    }
  }
}
