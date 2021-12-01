
#include <finite_field.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int64_t max(int64_t a, int64_t b) { return a > b ? a : b; }

static int64_t min(int64_t a, int64_t b) { return a < b ? a : b; }

static int64_t modulo(int64_t n, int64_t mod) { return (n % mod + mod) % mod; }

static void check_matching_prime_field(const struct prime_field_element *pfe1,
                                       const struct prime_field_element *pfe2) {
  if (pfe1->characteristic != pfe2->characteristic) {
    fprintf(stderr, "characteristic mismatch");
    exit(EXIT_FAILURE);
  }
}

struct prime_field_element *prime_field_element_new(int64_t capacity,
                                                    int64_t characteristic) {
  struct prime_field_element *pfe = malloc(sizeof(struct prime_field_element));
  if (pfe == NULL) {
    fprintf(stderr, "malloc");
    exit(EXIT_FAILURE);
  }

  pfe->coefficients = malloc(capacity * sizeof(int64_t));
  if (pfe->coefficients == NULL) {
    free(pfe);
    fprintf(stderr, "malloc");
    exit(EXIT_FAILURE);
  }

  pfe->capacity = capacity;
  pfe->characteristic = characteristic;
  pfe->degree = 0;
  memset(pfe->coefficients, 0, capacity * sizeof(int64_t));
  return pfe;
}

void prime_field_element_free(struct prime_field_element *pfe) {
  if (pfe != NULL) {
    free(pfe->coefficients);
    free(pfe);
  }
}

int64_t prime_field_element_degree(const struct prime_field_element *pfe) {
  for (int i = pfe->capacity - 1; i >= 0; i--) {
    if (pfe->coefficients[i] != 0) {
      return i;
    }
  }
  return 0; // TODO: negative degree for the zero polynomial
}

struct prime_field_element *
prime_field_element_from_array(int64_t *coefficients, int64_t capacity,
                               int64_t characteristic) {
  struct prime_field_element *pfe =
      prime_field_element_new(capacity, characteristic);

  for (int64_t i = 0; i < capacity; i++) {
    pfe->coefficients[i] = modulo(coefficients[i], pfe->characteristic);
  }

  pfe->degree = prime_field_element_degree(pfe);

  return pfe;
}

struct prime_field_element *
prime_field_element_copy(const struct prime_field_element *pfe) {
  return prime_field_element_from_array(pfe->coefficients, pfe->capacity,
                                        pfe->characteristic);
}

void prime_field_element_print(const struct prime_field_element *pfe) {
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

struct prime_field_element *
prime_field_element_add(const struct prime_field_element *pfe1,
                        const struct prime_field_element *pfe2) {
  check_matching_prime_field(pfe1, pfe2);

  struct prime_field_element *pfe = prime_field_element_new(
      max(pfe1->capacity, pfe2->capacity), pfe1->characteristic);

  int64_t i;
  for (i = 0; i < min(pfe1->degree, pfe2->degree); i++) {
    pfe->coefficients[i] = modulo(pfe1->coefficients[i] + pfe2->coefficients[i],
                                  pfe1->characteristic);
  }

  for (; i < max(pfe1->degree, pfe2->degree); i++) {
    pfe->coefficients[i] = (pfe1->degree > pfe2->degree)
                               ? pfe1->coefficients[i]
                               : pfe2->coefficients[i];
  }

  pfe->degree = prime_field_element_degree(pfe);

  return pfe;
}

struct prime_field_element *
prime_field_element_subtract(const struct prime_field_element *pfe1,
                             const struct prime_field_element *pfe2) {
  check_matching_prime_field(pfe1, pfe2);

  struct prime_field_element *pfe = prime_field_element_new(
      max(pfe1->capacity, pfe2->capacity), pfe1->characteristic);

  int64_t i;
  for (i = 0; i < min(pfe1->degree, pfe2->degree); i++) {
    pfe->coefficients[i] = modulo(pfe1->coefficients[i] - pfe2->coefficients[i],
                                  pfe1->characteristic);
  }

  for (; i < max(pfe1->degree, pfe2->degree); i++) {
    pfe->coefficients[i] = (pfe1->degree > pfe2->degree)
                               ? pfe1->coefficients[i]
                               : pfe2->coefficients[i];
  }

  pfe->degree = prime_field_element_degree(pfe);

  return pfe;
}

struct prime_field_element *
prime_field_element_multiplication(const struct prime_field_element *pfe1,
                                   const struct prime_field_element *pfe2) {
  check_matching_prime_field(pfe1, pfe2);

  struct prime_field_element *pfe = prime_field_element_new(
      pfe1->capacity * pfe2->capacity, pfe1->characteristic);

  for (int i = 0; i < pfe1->capacity; i++) {
    for (int j = 0; j < pfe2->capacity; j++) {
      pfe->coefficients[i + j] =
          modulo(pfe->coefficients[i + j] +
                     pfe1->coefficients[i] * pfe2->coefficients[j],
                 pfe1->characteristic);
    }
  }

  pfe->degree = prime_field_element_degree(pfe);

  return pfe;
}

static int64_t modular_inverse(int64_t a, int64_t b) {
  int64_t t, nt, r, nr, q, tmp;
  if (b < 0)
    b = -b;
  if (a < 0)
    a = b - (-a % b);

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
  if (r > 1)
    return -1;
  if (t < 0)
    t += b;
  return t;
}

void prime_field_element_remainder_mut(
    struct prime_field_element *pfe,
    const struct prime_field_element *divisor) {
  check_matching_prime_field(pfe, divisor);

  if (pfe->degree - divisor->degree < 0) {
    return;
  }

  struct prime_field_element *remainder = pfe;
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
  remainder->degree = prime_field_element_degree(remainder);
}

void prime_field_element_division(const struct prime_field_element *pfe1,
                                  const struct prime_field_element *pfe2,
                                  struct prime_field_element **remainder,
                                  struct prime_field_element **quotient) {
  check_matching_prime_field(pfe1, pfe2);

  int d = pfe1->degree - pfe2->degree;
  if (d < 0) {
    int64_t zero[1] = {0};
    *quotient = prime_field_element_from_array(zero, 1, pfe1->characteristic);
    *remainder = prime_field_element_copy(pfe1);
    return;
  }

  *quotient = prime_field_element_new(d + 1, pfe1->characteristic);
  *remainder = prime_field_element_copy(pfe1);

  int64_t i, j;
  int64_t ratio;

  for (i = pfe1->degree; i >= pfe2->degree; i--) {
    if ((*remainder)->coefficients[i] != 0) {
      ratio = (modular_inverse(pfe2->coefficients[pfe2->degree],
                               pfe1->characteristic) *
               (*remainder)->coefficients[i]);

      (*quotient)->coefficients[i - pfe2->degree] = ratio;
      (*remainder)->coefficients[i] = 0;

      for (j = 0; j < pfe2->degree; j++)
        (*remainder)->coefficients[i - pfe2->degree + j] =
            modulo((*remainder)->coefficients[i - pfe2->degree + j] -
                       (pfe2->coefficients[j] * ratio),
                   pfe1->characteristic);
    }
  }
  (*remainder)->degree = prime_field_element_degree(*remainder);
  (*quotient)->degree = prime_field_element_degree(*quotient);
}

struct prime_field_element *
prime_field_element_gcd(const struct prime_field_element *pfe1,
                        const struct prime_field_element *pfe2) {
  struct prime_field_element *r0, *r1, *t;

  if (pfe1->degree < pfe2->degree) {
    r0 = prime_field_element_copy(pfe2);
    r1 = prime_field_element_copy(pfe1);
  } else {
    r0 = prime_field_element_copy(pfe1);
    r1 = prime_field_element_copy(pfe2);
  }

  do {
    prime_field_element_remainder_mut(r0, r1);
    t = r0;
    r0 = r1;
    r1 = t;
  } while (!(r1->degree == 0 && *r1->coefficients == 0));

  free(r1);

  return r0;
}

struct prime_field_element *prime_field_element_gcd_extended(
    const struct prime_field_element *pfe1,
    const struct prime_field_element *pfe2,
    struct prime_field_element **bezout_coefficient_pfe1,
    struct prime_field_element **bezout_coefficient_pfe2) {
  check_matching_prime_field(pfe1, pfe2);

  struct prime_field_element *s = prime_field_element_new(
      pfe1->capacity + pfe2->capacity, pfe1->characteristic);
  struct prime_field_element *old_s = prime_field_element_new(
      pfe1->capacity + pfe2->capacity, pfe1->characteristic);
  old_s->coefficients[0] = 1;
  struct prime_field_element *t = prime_field_element_new(
      pfe1->capacity + pfe2->capacity, pfe1->characteristic);
  t->coefficients[0] = 1;
  struct prime_field_element *old_t = prime_field_element_new(
      pfe1->capacity + pfe2->capacity, pfe1->characteristic);
  struct prime_field_element *r = prime_field_element_copy(pfe2);
  struct prime_field_element *old_r = prime_field_element_copy(pfe1);

  struct prime_field_element *rem, *quot;

  while (!(r->degree == 0 && *r->coefficients == 0)) {
    prime_field_element_division(old_r, r, &rem, &quot);
    old_r = r;
    r = rem;
    old_s = s;
    s = prime_field_element_subtract(
        old_s, prime_field_element_multiplication(quot, s));
    old_t = t;
    t = prime_field_element_subtract(
        old_t, prime_field_element_multiplication(quot, t));
  }

  *bezout_coefficient_pfe1 = old_s;
  *bezout_coefficient_pfe2 = old_t;
  return old_r;
}

static void
check_matching_finite_field(const struct finite_field_element *ffe1,
                            const struct finite_field_element *ffe2) {
  check_matching_prime_field(ffe1->poly, ffe2->poly);
  check_matching_prime_field(ffe1->mod, ffe2->mod);
  if (memcmp(ffe1->mod, ffe2->mod, ffe1->mod->capacity + ffe2->mod->capacity)) {
    fprintf(stderr, "mismatching finite fields");
    exit(EXIT_FAILURE);
  }
}

struct finite_field_element *
finite_field_element_multiplication(const struct finite_field_element *ffe1,
                                    const struct finite_field_element *ffe2) {
  // check_matching_finite_field(ffe1, ffe2);

  struct finite_field_element *res =
      malloc(sizeof(struct finite_field_element));
  if (res == NULL) {
    fprintf(stderr, "malloc");
    exit(EXIT_FAILURE);
  }

  res->poly = prime_field_element_multiplication(ffe1->poly, ffe2->poly);
  prime_field_element_remainder_mut(res->poly, ffe1->mod);
  return res;
}