#include <fcntl.h>
#include <helpers.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

int64_t modulo(int64_t n, int64_t mod) { return (n % mod + mod) % mod; }

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

int64_t fast_pow(int64_t base, int64_t power) {
  int64_t res = 1;
  while (power > 0) {
    if ((power & 1) > 0) {
      res = res * base;
    }
    power >>= 1;
    base = base * base;
  }
  return res;
}

int64_t modular_inverse(int64_t a, int64_t b) {
  int64_t t, nt, r, nr, q, tmp;
  if (b < 0) {
    b = -b;
  }
  if (a < 0) {
    a = b - (-a % b);
  }
  t = 0;
  nt = 1;
  r = b;
  nr = a % b;

  while (nr != 0) {
    q = r / nr;

    tmp = nt;
    nt = t - q * nt;
    t = tmp;

    tmp = nr;
    nr = r - q * nr;
    r = tmp;
  }
  if (r > 1) {
    return -1;
  }
  if (t < 0) {
    t += b;
  }
  return t;
}

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

struct polynomial_in_Fp *
prime_field_element_power(const struct polynomial_in_Fp *pfe,
                          uint64_t power) {
  struct polynomial_in_Fp *res =
      polynomial_in_Fp_new(pfe->capacity, pfe->characteristic);
  res->coefficients[0] = 1;
  struct polynomial_in_Fp *pfe_copy = polynomial_in_Fp_copy(pfe);

  struct polynomial_in_Fp *tmp;

  while (power) {
    if (power & 1) {
      tmp = polynomial_in_Fp_multiplication(res, pfe_copy);
      polynomial_in_Fp_free(res);
      res = tmp;
    }
    tmp = polynomial_in_Fp_multiplication(pfe_copy, pfe_copy);
    polynomial_in_Fp_free(pfe_copy);
    pfe_copy = tmp;
    power >>= 1;
  }
  return res;
}

struct finite_field_element *
finite_field_element_power(const struct finite_field_element *ffe,
                           uint64_t power) {
  struct finite_field_element *res =
      malloc(sizeof(struct finite_field_element));
  res->poly =
      polynomial_in_Fp_new(ffe->poly->capacity, ffe->mod->characteristic);
  res->mod = polynomial_in_Fp_copy(ffe->mod);
  res->poly->coefficients[0] = 1;

  struct finite_field_element *ffe_copy = finite_field_element_copy(ffe);

  struct finite_field_element *tmp;
  
  while (power) {
    if (power & 1) {
      tmp = finite_field_element_multiplication(res, ffe_copy);
      finite_field_element_free(res);
      res = tmp;
    }
    tmp = finite_field_element_multiplication(ffe_copy, ffe_copy);
    finite_field_element_free(ffe_copy);
    ffe_copy = tmp;
    power >>= 1;
  }
  res->poly->degree = polynomial_in_Fp_degree(res->poly);
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