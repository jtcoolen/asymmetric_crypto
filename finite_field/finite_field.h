#include <stdint.h>

struct prime_field_element {
  int64_t *coefficients;
  int64_t capacity;
  int64_t characteristic;
  int64_t degree;
};

struct finite_field_element {
  struct prime_field_element *poly;
  struct prime_field_element *mod;
};

struct prime_field_element *prime_field_element_new(int64_t capacity,
                                                    int64_t characteristic);

void prime_field_element_element_free(struct prime_field_element *pfe);

int64_t prime_field_element_degree(const struct prime_field_element *pfe);

struct prime_field_element *
prime_field_element_from_array(int64_t *coefficients, int64_t len,
                               int64_t characteristic);

struct prime_field_element *
prime_field_element_copy(const struct prime_field_element *pfe);

void prime_field_element_print(const struct prime_field_element *pfe);

struct prime_field_element *
prime_field_element_add(const struct prime_field_element *pfe1,
                        const struct prime_field_element *pfe2);

struct prime_field_element *
prime_field_element_subtract(const struct prime_field_element *pfe1,
                             const struct prime_field_element *pfe2);

struct prime_field_element *
prime_field_element_multiplication(const struct prime_field_element *pfe1,
                                   const struct prime_field_element *pfe2);

void prime_field_element_remainder_mut(
    struct prime_field_element *pfe, const struct prime_field_element *divisor);

void prime_field_element_division(const struct prime_field_element *pfe1,
                                  const struct prime_field_element *pfe2,
                                  struct prime_field_element **remainder,
                                  struct prime_field_element **quotient);

struct prime_field_element *
prime_field_element_gcd(const struct prime_field_element *pfe1,
                        const struct prime_field_element *pfe2);

struct prime_field_element *prime_field_element_gcd_extended(
    const struct prime_field_element *pfe1,
    const struct prime_field_element *pfe2,
    struct prime_field_element **bezout_coefficient_pfe1,
    struct prime_field_element **bezout_coefficient_pfe2);

struct finite_field_element *
finite_field_element_multiplication(const struct finite_field_element *ffe1,
                                    const struct finite_field_element *ffe2);