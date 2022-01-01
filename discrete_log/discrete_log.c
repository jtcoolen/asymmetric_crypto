#include <blake2b.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <discrete_log.h>
#include <hashmap.h>
//#include <uthash.h>

struct babystep_pair {
  mpz_t exponent; // a
  mpz_t power;    // g^a
};

char *hash(const mpz_t a) {
  size_t s = mpz_sizeinbase(a, 10);
  char *str = malloc(s + 1);
  gmp_snprintf(str, s + 1, "%Zd", a);
  return str;
}

void babystep_giantstep(mpz_t cyclic_group_order, mpz_t generator, mpz_t h,
                        mpz_t log) {
  mpz_t group_order_sqrt, group_order_sqrt_neg, a, j, gen_pow, h_pow;
  mpz_inits(group_order_sqrt, group_order_sqrt_neg, a, j, gen_pow, h_pow, NULL);

  mpf_t sqrt;
  mpf_init(sqrt);
  mpf_set_z(sqrt, cyclic_group_order);
  mpf_sqrt(sqrt, sqrt);
  mpf_ceil(sqrt, sqrt);
  mpz_set_f(group_order_sqrt, sqrt);

  HASHMAP(char, struct babystep_pair) babystep_pairs;
  hashmap_init(&babystep_pairs, hashmap_hash_string, strcmp);

  mpz_set_ui(a, 1);

  struct babystep_pair *p;
  while (mpz_cmp(j, group_order_sqrt) <= 0) {
    p = malloc(sizeof *p);
    mpz_set(p->exponent, j);
    mpz_set(p->power, a);

    hashmap_put(&babystep_pairs, hash(a), p);

    mpz_mul(a, a, generator);
    mpz_mod(a, a, cyclic_group_order);
    mpz_add_ui(j, j, 1);
  }

  /*const mpz_t *k;

  hashmap_foreach(k, p, &babystep_pairs) {
    gmp_printf("\npairs: j= %Zd ; pow = %Zd \n", p->exponent, p->power);
  }*/

  mpz_neg(group_order_sqrt_neg, group_order_sqrt);
  mpz_powm(gen_pow, generator, group_order_sqrt_neg, cyclic_group_order);

  mpz_set(h_pow, h);
  mpz_set_ui(j, 0);
  while (mpz_cmp(j, group_order_sqrt) < 0) {
    p = hashmap_get(&babystep_pairs, hash(h_pow));
    if (p) {
      mpz_out_str(stdout, 10, p->exponent);
      mpz_mul(j, j, group_order_sqrt);
      mpz_add(log, j, p->exponent);
      hashmap_cleanup(&babystep_pairs);
      return;
    }
    mpz_mul(h_pow, h_pow, gen_pow);
    mpz_mod(h_pow, h_pow, cyclic_group_order);
    mpz_add_ui(j, j, 1);
  }
  hashmap_cleanup(&babystep_pairs);
}