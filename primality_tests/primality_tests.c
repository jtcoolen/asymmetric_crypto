#include <fcntl.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

int64_t modulo(int64_t n, int64_t mod) { return (n % mod + mod) % mod; }

int64_t randrange(int64_t lower, int64_t upper) {
  if (upper < lower) {
    fprintf(stderr, "wrong range");
    exit(EXIT_FAILURE);
  }
  int randomData = open("/dev/urandom", O_RDONLY);
  if (randomData < 0) {
    fprintf(stderr, "open");
    exit(EXIT_FAILURE);
  }
  char rd[8];
  int64_t rd_uint;
  ssize_t result = read(randomData, rd, sizeof rd);
  if (result < 0) {
    fprintf(stderr, "read");
    exit(EXIT_FAILURE);
  }
  memcpy(&rd_uint, rd, 8);
  return (modulo(rd_uint, (upper - lower + 1))) + lower;
}

int64_t pow_mod(int64_t base, int64_t power, int64_t mod) {
  int64_t res = 1;
  while (power > 0) {
    if ((power & 1) > 0) {
      res = modulo(res * base, mod);
    }
    power >>= 1;
    base = modulo(base * base, mod);
  }
  return res;
}

int64_t gcd(int64_t a, int64_t b) {
  int res;
  while ((a % b) > 0) {
    res = a % b;
    a = b;
    b = res;
  }
  return b;
}

int is_probably_prime_Fermat(int64_t n, size_t k) {
  int64_t a;
  for (size_t i = 0; i < k; i++) {
    a = randrange(2, n - 2);
    if (pow_mod(a, n - 1, n) != 1) {
      return 0;
    }
  }
  return 1;
}

#define SWAP(a, b) (((a) ^= (b)), ((b) ^= (a)), ((a) ^= (b)))

int Jacobi(int64_t a, int64_t n) {
  if (a >= n)
    a %= n;
  int result = 1;
  while (a) {
    while ((a & 1) == 0) {
      a >>= 1;
      if ((n & 7) == 3 || (n & 7) == 5)
        result = -result;
    }
    SWAP(a, n);
    if ((a & 3) == 3 && (n & 3) == 3)
      result = -result;
    a %= n;
  }
  if (n == 1)
    return result;
  return 0;
}

int is_probably_prime_Solovay_Strassen(int64_t n, size_t k) {
  int64_t a;
  for (size_t i = 0; i < k; i++) {
    a = randrange(2, n - 2);
    if (gcd(a, n) != 1) {
      return 0;
    }
    if (!((Jacobi(a, n) == -1 && pow_mod(a, (n - 1) >> 1, n) == n - 1) ||
          ((Jacobi(a, n) == 1 && pow_mod(a, (n - 1) >> 1, n) == 1)))) {
      return 0;
    }
  }
  return 1;
}

int main(void) {
  printf("is_probably_prime_Fermat(113)=%d\n",
         is_probably_prime_Fermat(113, 100));
  printf("is_probably_prime_Fermat(383)=%d\n",
         is_probably_prime_Fermat(383, 100));
  printf("is_probably_prime_Fermat(187)=%d\n",
         is_probably_prime_Fermat(187, 100));
  printf("is_probably_prime_Fermat(91)=%d\n",
         is_probably_prime_Fermat(91, 100));

  printf("is_probably_prime_Solovay_Strassen(113)=%d\n",
         is_probably_prime_Solovay_Strassen(113, 100));
  printf("is_probably_prime_Solovay_Strassen(383)=%d\n",
         is_probably_prime_Solovay_Strassen(383, 100));
  printf("is_probably_prime_Solovay_Strassen(187)=%d\n",
         is_probably_prime_Solovay_Strassen(187, 100));
  printf("is_probably_prime_Solovay_Strassen(91)=%d\n",
         is_probably_prime_Solovay_Strassen(91, 100));
}
