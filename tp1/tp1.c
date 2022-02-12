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

int main(void) {
  printf("%ld", int_sqrt(50));
  return 0;
}
