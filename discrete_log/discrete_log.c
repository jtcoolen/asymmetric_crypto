#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <discrete_log.h>
#include <hashmap.h>

struct babystep_pair {
  mpz_t exponent; // a
  mpz_t power;    // g^a
};

void babystep_giantstep(mpz_t cyclic_group_order, mpz_t generator, mpz_t h,
                        mpz_t log) {
  mpz_t group_order_sqrt, group_order_sqrt_neg, a, j, gen_pow, h_pow;
  mpz_inits(group_order_sqrt, group_order_sqrt_neg, a, j, gen_pow, h_pow, NULL);

  mpz_sqrt(group_order_sqrt, cyclic_group_order);

  HASHMAP(char, struct babystep_pair) babystep_pairs;

  hashmap_init(&babystep_pairs, hashmap_hash_string, strcmp);

  mpz_set_ui(a, 1);
  char key[64];
  struct babystep_pair *p;
  while (mpz_cmp(j, group_order_sqrt) < 0) {
    memset(key, 0, 64);
    printf("len=%d", gmp_snprintf(key, 64, "%Za", a));
    p = malloc(sizeof *p);
    mpz_set(p->exponent, j);
    mpz_set(p->power, a);
    hashmap_put(&babystep_pairs, key, p);
    mpz_mul(a, a, generator);
    mpz_mod(a, a, cyclic_group_order);
    mpz_add_ui(j, j, 1);
  }

  const char *k;

  hashmap_foreach(k, p, &babystep_pairs) {
    gmp_printf("\npairs[%s]: j= %Zd ; pow = %Zd \n", key, p->exponent, p->power);
  }

  mpz_neg(group_order_sqrt_neg, group_order_sqrt);
  mpz_powm(gen_pow, a, group_order_sqrt_neg, cyclic_group_order);

  mpz_set(h_pow, h);
  mpz_set_ui(j, 0);
  while (mpz_cmp(j, group_order_sqrt) < 0) {
    memset(key, 0, 64);
    gmp_printf("%Za \n", h_pow);
    printf("len=%d", gmp_snprintf(key, 64, "%Za", h_pow));
    p = hashmap_get(&babystep_pairs, key);
    if (p) {
      printf("\n");
      mpz_out_str(stdout, 10, p->power);
      mpz_mul(j, j, group_order_sqrt);
      mpz_add(log, j, p->exponent);
      return;
    }
    mpz_mul(h_pow, h_pow, gen_pow);
    mpz_add_ui(j, j, 1);
  }
}