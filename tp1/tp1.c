#include <stdint.h>
#include <stdio.h>

int log2(uint64_t n) {
  int d = 0;
  while (n) {
    d++;
    n >>= 1;
  }
  return d;
}

uint64_t int_sqrt(uint64_t n) {
  uint64_t sqrt = 0;
  uint64_t approx;
  for (int i = (log2(n)>>1) + 1; i >= 0; i--) {
    approx = sqrt ^ (1 << i);
    if (approx * approx <= n) {
      sqrt ^= (1 << i);
    }
  }
  return sqrt;
}

int64_t abs(int64_t n) { return (n < 0) ? -n : n; }

int64_t gcd(int64_t a, int64_t b) {
  int64_t t;
  while (b != 0) {
    t = b;
    // On choisit le plus petit reste
    // On gagne un facteur 2 car à chaque itération le reste choisi r_i est < |r_(i-1)|/2
    // en valeur absolue
    if (abs(a % b) > abs(a % b - b)) {
      b = a % b - b;
    } else {
      b = a % b;
    }
    a = t;
  }
  return a;
}

int main(void) {
  printf("int_sqrt(50)=%ld\n", int_sqrt(50));
  printf("gcd(50, 125)=%ld", gcd(50, 125));
  return 0;
}
