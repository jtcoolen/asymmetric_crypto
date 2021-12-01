#include <stdint.h>

uint64_t pow_mod(uint64_t base, uint64_t power, uint64_t mod) {
  uint64_t res = 1;
  while (power > 0) {
    if ((power & 1) > 0) {
      res = (res * base) % mod;
    }
    power >>= 1;
    base = (base * base) % mod;
  }
  return res;
}

int64_t find_factor(int64_t a) {
  
}

int legendre_symbol(int64_t a, int64_t p) {  
  
  if (a > p) {
    a %= p;
  }

  if (a == 1) {
    return 1;
  }

  if (a == 2) {
    if (p % 8 == 1 || p % 8 == 7) {
      return 1;
    }
    if (p % 8 == 3 || p % 8 == 5) {
      return -1;
    }
  }

  if (!is_prime(a)) {
    return 
  }

  // a is prime
  // quadratic reciprocity law
  if (a % 4 == 3 && p % 4 == 3) {
    return -1 * legendre_symbol(p, a);
  }

  return legendre_symbol(p, a);

}