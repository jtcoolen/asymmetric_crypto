#include <stdint.h>
#include <stdio.h>

uint64_t log(uint64_t n) {
  uint64_t d = 0;
  while (n) {
    d++;
    n >>= 1;
  }
  return d;
}

uint64_t int_sqrt(uint64_t n) {
  uint64_t sqrt = (1 << log(n));
  uint64_t approx;
  for (uint64_t i = log(n) + 1; i >= 0; i--) {
    approx = sqrt * sqrt;
    if (n == approx) {
      return sqrt;
    }
    if (n < approx) {
      sqrt = sqrt & ~(1ull << i);
      sqrt = sqrt & ~(1ull << (i - 1));
      sqrt |= (1ull << (i - 1));
    } else {
      sqrt = sqrt & ~(1ull << (i - 1));
      sqrt |= (1ull << (i - 1));
    }
  }
  return 0;
}

int main(void) {
  printf("%ld", int_sqrt(50));
  return 0;
}