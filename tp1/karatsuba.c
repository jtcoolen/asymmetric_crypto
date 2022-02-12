#include <inttypes.h>
#include <stdlib.h>
#include <stdio.h>

/*int
size_base2(uint64_t n)
{
  int s = 0;
  while (n) {
    s++;
    n >>= 1;
  }
  return s;
  }*/
int
size_base2(uint64_t n, int r)
{
  if (!n) return r;
  n >>= 1;
  return size_base2(n, ++r);
}

uint64_t
min(uint64_t n, uint64_t m)
{
  return (n < m) ? n : m;
}

uint64_t
karatsuba(uint64_t n, uint64_t m)
{
  if ((n < 2) || (m < 2))
    return n * m;

  int s = min(size_base2(n, 0), size_base2(m, 0));

  s = (s / 2);

  uint64_t h1 = (n >> s);
  uint64_t l1 = n ^ (h1 << s);

  uint64_t h2 = (m >> s);
  uint64_t l2 = m ^ (h2 << s);

  uint64_t z0 = karatsuba(l1, l2);
  uint64_t z1 = karatsuba(h1 - l1, h2 - l2);
  uint64_t z2 = karatsuba(h1, h2);

  return z2 * (1 << (2 * s)) + (z2 + z0 - z1) * (1 << s) + z0;
}

int
main(void)
{
  uint64_t n, m;
  //for(int i=0; i< 10000000;i++) {
  n = 932023;
  m = 105978;
  //printf("n=%ld, m=%ld\n", n, m);
  uint64_t r = karatsuba(n, m);
  printf("\nr=%ld", r);
  //}
  return 0;
}
