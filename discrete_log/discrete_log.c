#include <discrete_log.h>
#include <hashmap.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

struct babystep_pair {
  mpz_t exponent; // a
  mpz_t power;    // g^a
};

void babystep_giantstep(mpz_t cyclic_group_order, mpz_t generator, mpz_t h) {
  mpz_t group_order_sqrt, a, j;
  mpz_inits(group_order_sqrt, a, j, NULL);

  mpz_sqrt(group_order_sqrt, cyclic_group_order);

  HASHMAP(char, struct babystep_pair) babystep_pairs;

  hashmap_init(&babystep_pairs, hashmap_hash_string, strcmp);

  mpz_set_ui(a, 1);
  char key[64];
  struct babystep_pair *p;
  while (mpz_cmp(j, group_order_sqrt) < 0) {
    printf("len=%d", gmp_snprintf(key, 64, "%Za", a));
    p = malloc(sizeof *p);
    mpz_set(p->exponent, j);
    mpz_set(p->power, a);
    hashmap_put(&babystep_pairs, key, p);
    mpz_mul(a, a, generator);
    mpz_mod(a, a, cyclic_group_order);
    mpz_add_ui(j, j, 1);
  }
}