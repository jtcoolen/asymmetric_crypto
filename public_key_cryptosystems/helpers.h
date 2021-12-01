#ifndef HELPERS_H
#define HELPERS_H

#include <finite_field.h>
#include <stdint.h>

int64_t modulo(int64_t n, int64_t mod);

int64_t pow_mod(int64_t base, int64_t power, int64_t mod);

int64_t fast_pow(int64_t base, int64_t power);

int64_t modular_inverse(int64_t a, int64_t b);

int64_t randrange(int64_t lower, int64_t upper);

struct polynomial_in_Fp *
prime_field_element_power(const struct polynomial_in_Fp *pfe,
                          uint64_t power);

struct finite_field_element *
finite_field_element_power(const struct finite_field_element *ffe,
                           uint64_t power);

int64_t gcd(int64_t a, int64_t b);

#endif