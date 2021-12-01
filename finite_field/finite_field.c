#include <finite_field.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int64_t max(int64_t a, int64_t b) { return a > b ? a : b; }

static int64_t min(int64_t a, int64_t b) { return a < b ? a : b; }

static int64_t modulo(int64_t n, int64_t mod) { return (n % mod + mod) % mod; }

static uint64_t exponentiation(uint64_t base, uint64_t power) {
  uint64_t res = 1;
  while (power) {
    if (power & 1) {
      res *= base;
    }
    base *= base;
    power >>= 1;
  }
  return res;
}

static void check_matching_polynomial_in_Fp(const struct polynomial_in_Fp *pfe1,
                                       const struct polynomial_in_Fp *pfe2) {
  if (pfe1->characteristic != pfe2->characteristic) {
    fprintf(stderr, "characteristic mismatch");
    exit(EXIT_FAILURE);
  }
}

struct polynomial_in_Fp *polynomial_in_Fp_new(uint64_t capacity,
                                                    uint64_t characteristic) {
  struct polynomial_in_Fp *pfe = malloc(sizeof(struct polynomial_in_Fp));
  if (pfe == NULL) {
    fprintf(stderr, "malloc new 1 ");
    exit(EXIT_FAILURE);
  }

  pfe->coefficients = malloc(capacity * sizeof(int64_t));
  if (pfe->coefficients == NULL) {
    free(pfe);
    fprintf(stderr, "malloc new 2 capa = %ld", capacity);
    exit(EXIT_FAILURE);
  }

  pfe->capacity = capacity;
  pfe->characteristic = characteristic;
  pfe->degree = 0;
  memset(pfe->coefficients, 0, capacity * sizeof(int64_t));
  return pfe;
}

void polynomial_in_Fp_free(struct polynomial_in_Fp *pfe) {
  if (pfe != NULL) {
    free(pfe->coefficients);
    free(pfe);
  }
}

int64_t polynomial_in_Fp_degree(const struct polynomial_in_Fp *pfe) {
  for (int i = pfe->capacity - 1; i >= 0; i--) {
    if (pfe->coefficients[i] != 0) {
      return i;
    }
  }
  return 0; // TODO: negative degree for the zero polynomial
}

struct polynomial_in_Fp *
polynomial_in_Fp_from_array(int64_t *coefficients, uint64_t capacity,
                               uint64_t characteristic) {
  struct polynomial_in_Fp *pfe =
      polynomial_in_Fp_new(capacity, characteristic);

  for (uint64_t i = 0; i < capacity; i++) {
    pfe->coefficients[i] = modulo(coefficients[i], pfe->characteristic);
  }

  pfe->degree = polynomial_in_Fp_degree(pfe);

  return pfe;
}

struct polynomial_in_Fp *
polynomial_in_Fp_copy(const struct polynomial_in_Fp *pfe) {
  return polynomial_in_Fp_from_array(pfe->coefficients, pfe->capacity,
                                        pfe->characteristic);
}

void polynomial_in_Fp_print(const struct polynomial_in_Fp *pfe) {
  if (pfe->degree == 0) {
    printf("%ld", pfe->coefficients[0]);
    return;
  }

  if (pfe->coefficients[0] != 0) {
    printf("%ld+", pfe->coefficients[0]);
  }

  for (int64_t i = 1; i < pfe->degree; i++) {
    if (pfe->coefficients[i] != 0)
      printf("%ld*X^%ld+", pfe->coefficients[i], i);
  }

  if (pfe->coefficients[pfe->degree] != 0) {
    printf("%ld*X^%ld", pfe->coefficients[pfe->degree], pfe->degree);
  }
}

struct polynomial_in_Fp *
polynomial_in_Fp_add(const struct polynomial_in_Fp *pfe1,
                        const struct polynomial_in_Fp *pfe2) {
  check_matching_polynomial_in_Fp(pfe1, pfe2);

  struct polynomial_in_Fp *pfe = polynomial_in_Fp_new(
      max(pfe1->capacity, pfe2->capacity), pfe1->characteristic);

  int64_t i;
  for (i = 0; i <= min(pfe1->degree, pfe2->degree); i++) {
    pfe->coefficients[i] = modulo(pfe1->coefficients[i] + pfe2->coefficients[i],
                                  pfe1->characteristic);
  }

  for (; i <= max(pfe1->degree, pfe2->degree); i++) {
    pfe->coefficients[i] = (pfe1->degree > pfe2->degree)
                               ? pfe1->coefficients[i]
                               : pfe2->coefficients[i];
  }

  pfe->degree = polynomial_in_Fp_degree(pfe);

  return pfe;
}

struct polynomial_in_Fp *
polynomial_in_Fp_subtract(const struct polynomial_in_Fp *pfe1,
                             const struct polynomial_in_Fp *pfe2) {
  check_matching_polynomial_in_Fp(pfe1, pfe2);

  struct polynomial_in_Fp *pfe = polynomial_in_Fp_new(
      max(pfe1->capacity, pfe2->capacity), pfe1->characteristic);

  int64_t i;
  for (i = 0; i <= min(pfe1->degree, pfe2->degree); i++) {
    pfe->coefficients[i] = modulo(pfe1->coefficients[i] - pfe2->coefficients[i],
                                  pfe1->characteristic);
  }

  for (; i <= max(pfe1->degree, pfe2->degree); i++) {
    pfe->coefficients[i] =
        (pfe1->degree > pfe2->degree)
            ? pfe1->coefficients[i]
            : modulo(-pfe2->coefficients[i], pfe1->characteristic);
  }

  pfe->degree = polynomial_in_Fp_degree(pfe);

  return pfe;
}

struct polynomial_in_Fp *
polynomial_in_Fp_multiplication(const struct polynomial_in_Fp *pfe1,
                                   const struct polynomial_in_Fp *pfe2) {
  check_matching_polynomial_in_Fp(pfe1, pfe2);

  int64_t pfe1_deg = polynomial_in_Fp_degree(pfe1);
  int64_t pfe2_deg = polynomial_in_Fp_degree(pfe2);

  struct polynomial_in_Fp *pfe = polynomial_in_Fp_new(
    (pfe1_deg + 50) * (pfe2_deg + 50), pfe1->characteristic);

  for (int64_t i = 0; i <= pfe1_deg; i++) {
    for (int64_t j = 0; j <= pfe2_deg; j++) {
      pfe->coefficients[i + j] =
          modulo(pfe->coefficients[i + j] +
                     pfe1->coefficients[i] * pfe2->coefficients[j],
                 pfe1->characteristic);
    }
  }

  pfe->degree = polynomial_in_Fp_degree(pfe);

  return pfe;
}

static int64_t modular_inverse(int64_t a, int64_t b) {
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

void prime_field_element_remainder(struct polynomial_in_Fp *pfe,
                                   const struct polynomial_in_Fp *divisor) {
  check_matching_polynomial_in_Fp(pfe, divisor);

  if (pfe->degree - divisor->degree < 0) {
    return;
  }

  struct polynomial_in_Fp *remainder = pfe;
  int64_t i, j;
  int64_t ratio;

  for (i = pfe->degree; i >= divisor->degree; i--) {
    if (remainder->coefficients[i] != 0) {
      ratio = modular_inverse(divisor->coefficients[divisor->degree],
                              pfe->characteristic) *
              remainder->coefficients[i];
      remainder->coefficients[i] = 0;

      for (j = 0; j < divisor->degree; j++)
        remainder->coefficients[i - divisor->degree + j] =
            modulo(remainder->coefficients[i - divisor->degree + j] -
                       divisor->coefficients[j] * ratio,
                   pfe->characteristic);
    }
  }
  remainder->degree = polynomial_in_Fp_degree(remainder);
}

void polynomial_in_Fp_division(const struct polynomial_in_Fp *pfe1,
                                  const struct polynomial_in_Fp *pfe2,
                                  struct polynomial_in_Fp **remainder,
                                  struct polynomial_in_Fp **quotient) {
  check_matching_polynomial_in_Fp(pfe1, pfe2);

  int64_t d = pfe1->degree - pfe2->degree;
  if (d < 0) {
    int64_t zero[1] = {0};
    *quotient = polynomial_in_Fp_from_array(zero, 1, pfe1->characteristic);
    *remainder = polynomial_in_Fp_copy(pfe1);
    return;
  }

  *quotient = polynomial_in_Fp_new(d + 1, pfe1->characteristic);
  *remainder = polynomial_in_Fp_copy(pfe1);

  int64_t i, j;
  int64_t ratio;

  for (i = pfe1->degree; i >= pfe2->degree; i--) {
    if ((*remainder)->coefficients[i] != 0) {
      ratio = modular_inverse(pfe2->coefficients[pfe2->degree],
                              pfe1->characteristic) *
              (*remainder)->coefficients[i];

      (*quotient)->coefficients[i - pfe2->degree] = ratio;
      (*remainder)->coefficients[i] = 0;

      for (j = 0; j < pfe2->degree; j++)
        (*remainder)->coefficients[i - pfe2->degree + j] =
            modulo((*remainder)->coefficients[i - pfe2->degree + j] -
                       (pfe2->coefficients[j] * ratio),
                   pfe1->characteristic);
    }
  }
  (*remainder)->degree = polynomial_in_Fp_degree(*remainder);
  (*quotient)->degree = polynomial_in_Fp_degree(*quotient);
}

struct polynomial_in_Fp *
polynomial_in_Fp_gcd(const struct polynomial_in_Fp *pfe1,
                        const struct polynomial_in_Fp *pfe2) {
  struct polynomial_in_Fp *r0, *r1, *t;

  if (pfe1->degree < pfe2->degree) {
    r0 = polynomial_in_Fp_copy(pfe2);
    r1 = polynomial_in_Fp_copy(pfe1);
  } else {
    r0 = polynomial_in_Fp_copy(pfe1);
    r1 = polynomial_in_Fp_copy(pfe2);
  }

  do {
    prime_field_element_remainder(r0, r1);
    t = r0;
    r0 = r1;
    r1 = t;
  } while (!(r1->degree == 0 && *r1->coefficients == 0));

  polynomial_in_Fp_free(r1);

  return r0;
}

struct polynomial_in_Fp *prime_field_element_gcd_extended(
    const struct polynomial_in_Fp *pfe1,
    const struct polynomial_in_Fp *pfe2,
    struct polynomial_in_Fp **bezout_coefficient_pfe1,
    struct polynomial_in_Fp **bezout_coefficient_pfe2) {
  check_matching_polynomial_in_Fp(pfe1, pfe2);

  struct polynomial_in_Fp *s = polynomial_in_Fp_new(
      max(pfe1->capacity, pfe2->capacity), pfe1->characteristic);

  struct polynomial_in_Fp *old_s = polynomial_in_Fp_new(
      max(pfe1->capacity, pfe2->capacity), pfe1->characteristic);
  old_s->coefficients[0] = 1;

  struct polynomial_in_Fp *t = polynomial_in_Fp_new(
      max(pfe1->capacity, pfe2->capacity), pfe1->characteristic);
  t->coefficients[0] = 1;

  struct polynomial_in_Fp *old_t = polynomial_in_Fp_new(
      max(pfe1->capacity, pfe2->capacity), pfe1->characteristic);

  struct polynomial_in_Fp *r = polynomial_in_Fp_copy(pfe2);

  struct polynomial_in_Fp *old_r = polynomial_in_Fp_copy(pfe1);

  struct polynomial_in_Fp *rem, *quot, *mul, *old_s_tmp, *old_t_tmp;

  while (!(r->degree == 0 && *r->coefficients == 0)) {
    polynomial_in_Fp_division(old_r, r, &rem, &quot);

    polynomial_in_Fp_free(old_r);
    old_r = r;

    r = rem;

    old_s_tmp = old_s;
    old_s = s;

    mul = polynomial_in_Fp_multiplication(quot, s);
    s = polynomial_in_Fp_subtract(old_s_tmp, mul);
    polynomial_in_Fp_free(mul);
    polynomial_in_Fp_free(old_s_tmp);

    old_t_tmp = old_t;
    old_t = t;

    mul = polynomial_in_Fp_multiplication(quot, t);
    t = polynomial_in_Fp_subtract(old_t_tmp, mul);
    polynomial_in_Fp_free(mul);
    polynomial_in_Fp_free(old_t_tmp);

    polynomial_in_Fp_free(quot);
  }
  polynomial_in_Fp_free(rem);
  polynomial_in_Fp_free(t);
  polynomial_in_Fp_free(s);

  *bezout_coefficient_pfe1 = old_s;
  *bezout_coefficient_pfe2 = old_t;
  return old_r;
}

static void
check_matching_finite_field(const struct finite_field_element *ffe1,
                            const struct finite_field_element *ffe2) {
  check_matching_polynomial_in_Fp(ffe1->poly, ffe2->poly);
  check_matching_polynomial_in_Fp(ffe1->mod, ffe2->mod);
  if (memcmp(ffe1->mod->coefficients, ffe2->mod->coefficients,
             min(ffe1->mod->capacity, ffe2->mod->capacity))) {
    fprintf(stderr, "mismatching finite fields");
    exit(EXIT_FAILURE);
  }
}

void finite_field_element_free(struct finite_field_element *ffe) {
  polynomial_in_Fp_free(ffe->poly);
  polynomial_in_Fp_free(ffe->mod);
  free(ffe);
}

struct finite_field_element *
finite_field_element_new(const struct polynomial_in_Fp *poly,
                         const struct polynomial_in_Fp *mod) {
  struct finite_field_element *ffe =
      malloc(sizeof(struct finite_field_element));
  ffe->poly = polynomial_in_Fp_copy(poly);
  ffe->mod = polynomial_in_Fp_copy(mod);
  ffe->multiplicative_group_order =
      exponentiation(mod->characteristic, mod->degree) - 1;
  return ffe;
}

struct finite_field_element *
finite_field_element_copy(const struct finite_field_element *ffe) {
  struct finite_field_element *ffe_copy =
      malloc(sizeof(struct finite_field_element));
  ffe_copy->mod = polynomial_in_Fp_copy(ffe->mod);
  ffe_copy->poly = polynomial_in_Fp_copy(ffe->poly);
  ffe_copy->multiplicative_group_order = ffe->multiplicative_group_order;
  return ffe_copy;
}

struct finite_field_element *
finite_field_element_multiplication(const struct finite_field_element *ffe1,
                                    const struct finite_field_element *ffe2) {
  check_matching_finite_field(ffe1, ffe2);

  struct finite_field_element *res =
      malloc(sizeof(struct finite_field_element));
  if (res == NULL) {
    fprintf(stderr, "malloc2");
    exit(EXIT_FAILURE);
  }

  res->poly = polynomial_in_Fp_multiplication(ffe1->poly, ffe2->poly);
  res->mod = polynomial_in_Fp_copy(ffe1->mod);
  prime_field_element_remainder(res->poly, ffe1->mod);
  return res;
}

struct polynomial_in_Fp *prime_field_element_scalar_multiplication(
    int64_t scalar, const struct polynomial_in_Fp *pfe) {
  struct polynomial_in_Fp *prod = polynomial_in_Fp_copy(pfe);
  for (int64_t i = 0; i <= pfe->degree; i++) {
    prod->coefficients[i] =
        modulo(prod->coefficients[i] * scalar, prod->characteristic);
  }
  return prod;
}

struct finite_field_element *
finite_field_element_inverse(const struct finite_field_element *ffe) {
  struct finite_field_element *inverse =
      malloc(sizeof(struct finite_field_element));
  if (inverse == NULL) {
    fprintf(stderr, "malloc");
    exit(EXIT_FAILURE);
  }
  // ffe->poly*a + ffe->mod b = gcd => ffe->poly*a = gcd [mod ffe->mod],
  // gcd invertible (constant) since ffe->mod is irreducible
  // ffe->poly*a*gcd^{-1} = 1 [mod ffe->mod]
  // thus inverse->poly = a*gcd^{-1}
  struct polynomial_in_Fp *a, *b, *gcd;
  gcd = prime_field_element_gcd_extended(ffe->poly, ffe->mod, &a, &b);

  inverse->poly = prime_field_element_scalar_multiplication(
      modular_inverse(gcd->coefficients[0], ffe->mod->characteristic), a);

  polynomial_in_Fp_free(gcd);
  polynomial_in_Fp_free(a);
  polynomial_in_Fp_free(b);

  inverse->mod = polynomial_in_Fp_copy(ffe->mod);

  return inverse;
}

struct finite_field_element *
finite_field_element_division(const struct finite_field_element *ffe1,
                              const struct finite_field_element *ffe2) {
  struct finite_field_element *ffe2_inv = finite_field_element_inverse(ffe2);
  struct finite_field_element *div =
      finite_field_element_multiplication(ffe1, ffe2_inv);
  finite_field_element_free(ffe2_inv);
  return div;
}

int finite_field_element_is_generator(const struct finite_field_element *ffe) {
  struct finite_field_element *pow = finite_field_element_new(ffe->poly, ffe->mod);
  uint64_t i = 1;
  while (!(pow->poly->degree != 0 && *pow->poly->coefficients == 0)) {
    struct finite_field_element *tmp = finite_field_element_multiplication(pow, ffe);
    finite_field_element_free(pow);
    pow = tmp;
    i++;
  }
  return i == ffe->multiplicative_group_order; 
}