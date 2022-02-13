#include <blake2b.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <discrete_log.h>
#include <hashmap.h>

struct babystep_pair {
  char *key;
  mpz_t exponent; // a
  mpz_t power;    // g^a
};


// Clé d'une entrée dans la table associative, définie comme la représentation en base 10 de l'entier
char *get_key(const mpz_t a) {
  size_t s = mpz_sizeinbase(a, 10);
  char *str = malloc(s + 1);
  gmp_snprintf(str, s + 1, "%Zd", a);
  return str;
}


// Version O(sqrt(N)) avec N l'ordre du groupe cyclique
// (O(n/t) avec t=sqrt(N))
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
  mpf_clear(sqrt);

  HASHMAP(char, struct babystep_pair) babystep_pairs;
  hashmap_init(&babystep_pairs, hashmap_hash_string, strcmp);

  mpz_set_ui(a, 1);

  // Précalcul O(t)=O(sqrt(N))
  // On calcule les paires (j, g^j)
  struct babystep_pair *p;
  while (mpz_cmp(j, group_order_sqrt) <= 0) {
    p = malloc(sizeof *p);
    if (p == NULL) {
      return;
    }
    mpz_inits(p->exponent, p->power, NULL);
    mpz_set(p->exponent, j);
    mpz_set(p->power, a);
    p->key = get_key(a);

    hashmap_put(&babystep_pairs, p->key, p);

    mpz_mul(a, a, generator);
    mpz_mod(a, a, cyclic_group_order);
    mpz_add_ui(j, j, 1);
  }

  // Affichage débogage
  /*const mpz_t *k;

  hashmap_foreach(k, p, &babystep_pairs) {
    gmp_printf("\npairs: j= %Zd ; pow = %Zd \n", p->exponent, p->power);
  }*/

  mpz_neg(group_order_sqrt_neg, group_order_sqrt);
  mpz_powm(gen_pow, generator, group_order_sqrt_neg, cyclic_group_order);

  // Calcul du log discret de h en base g (=variable 'generator')
  mpz_set(h_pow, h);
  mpz_set_ui(j, 0);
  char *pow_key;
  // On calcule h*g^(-j*t)=h*g^(-t)*g^(-t)* ... * g^(-t) (j <= sqrt(N)-1 fois)
  while (mpz_cmp(j, group_order_sqrt) < 0) {
    pow_key = get_key(h_pow);
    p = hashmap_get(&babystep_pairs, pow_key);
    free(pow_key);
    if (p) { // Collision h g^(-t*j) = p (=g^k) => h = g^(j+t*k) avec t=sqrt(N)
      mpz_mul(j, j, group_order_sqrt);
      // on a donc x = j * sqrt(N) + k, (k est l'exposant de p=g^k)
      mpz_add(log, j, p->exponent);

      // Suppression de la table de hachage
      hashmap_foreach_data(p, &babystep_pairs) {
        free(p->key);
        free(p);
      }
      hashmap_cleanup(&babystep_pairs);
      mpz_clears(group_order_sqrt, group_order_sqrt_neg, a, j, gen_pow, h_pow,
                 NULL);

      return;
    }
    // h=h*(g^(-t))
    mpz_mul(h_pow, h_pow, gen_pow);
    mpz_mod(h_pow, h_pow, cyclic_group_order);
    mpz_add_ui(j, j, 1);
  }

  // Suppression de la table de hachage
  hashmap_foreach_data(p, &babystep_pairs) {
    free(p->key);
    free(p);
  }
  hashmap_cleanup(&babystep_pairs);
  mpz_clears(group_order_sqrt, group_order_sqrt_neg, a, j, gen_pow, h_pow,
             NULL);
}

// Résout la congruence a * x = b [mod] pour x
int solve_congruence(mpz_t a, mpz_t b, mpz_t mod, mpz_t res) {
  mpz_t gcd, m, s, t, zero;
  mpz_inits(gcd, m, s, t, zero, NULL);
  mpz_gcdext(gcd, s, t, a, mod);

  mpz_mod(m, b, gcd);
  // une solution existe ssi b est divisible par le pgcd de a et mod
  if (mpz_cmp(m, zero) != 0) {
    return -1;
  }

  // les solutions sont de la forme: for(i=1,gcd-1, print((b/gcd)*s+i*(mod/gcd)) )
  mpz_divexact(res, b, gcd);
  mpz_mul(res, res, s);
  mpz_mod(res, res, mod);
  return 1;
}

unsigned long partition(mpz_t y) {
  mpz_t res;
  mpz_init(res);
  mpz_mod_ui(res, y, 3);
  unsigned long res_ui = mpz_get_ui(res);
  mpz_clear(res);
  return res_ui;
}

// fonction pseudoaléatoire pour la marche aléatoire dans un groupe cyclique d'ordre cyclic_group_order
void f_seq_terms(mpz_t y, mpz_t a, mpz_t b, mpz_t cyclic_group_order,
                 mpz_t generator, mpz_t pow) {
  mpz_t ord;
  mpz_init(ord);
  mpz_sub_ui(ord, cyclic_group_order, 1);
  switch (partition(y)) { // le groupe ambiant est partitionné en 3 sous-ensembles, on récupère la classe de l'élément y
  case 0:
    mpz_mul(y, y, y);
    mpz_mod(y, y, cyclic_group_order);
    mpz_mul_ui(a, a, 2);
    mpz_mod(a, a, ord);
    mpz_mul_ui(b, b, 2);
    mpz_mod(b, b, ord);
    break;
  case 1:
    mpz_mul(y, y, generator);
    mpz_mod(y, y, cyclic_group_order);
    mpz_add_ui(a, a, 1);
    mpz_mod(a, a, ord);
    break;
  case 2:
    mpz_mul(y, y, pow);
    mpz_mod(y, y, cyclic_group_order);
    mpz_add_ui(b, b, 1);
    mpz_mod(b, b, ord);
    break;
  }
}

// Rho-Pollard avec mémoire constante
// la complexité atteint la borne de Shoup et est donc optimale (algo non déterministe)
void pollard_rho(mpz_t cyclic_group_order, mpz_t generator, mpz_t h,
                 mpz_t log) {

  mpz_t yk, ak, bk, y2k, a2k, b2k;
  mpz_inits(yk, ak, bk, y2k, a2k, b2k, NULL);
  mpz_set_ui(yk, 1);
  mpz_set_ui(ak, 0);
  mpz_set_ui(bk, 0);

  mpz_set(y2k, yk);
  mpz_set(a2k, ak);
  mpz_set(b2k, bk);

  while (1) {
    // pas de tortue
    f_seq_terms(yk, ak, bk, cyclic_group_order, generator, h);

    // pas de lièvre
    f_seq_terms(y2k, a2k, b2k, cyclic_group_order, generator, h);
    f_seq_terms(y2k, a2k, b2k, cyclic_group_order, generator, h);

    // recherche d'une collision y_k = y_{2k}
    if (mpz_cmp(yk, y2k) == 0) {
      // on obtient une collision  g^(a2k-ak) = h^(bk-b2k)
      // dont on tire : h=g^x avec (b2k-bk)x = (ak-a2k) mod (n-1)
      mpz_sub(ak, ak, a2k);
      mpz_sub(bk, b2k, bk);
      mpz_t mod;
      mpz_init(mod);
      mpz_sub_ui(mod, cyclic_group_order, 1);
      solve_congruence(bk, ak, mod, log); // On résout la congruence pour x
      return;
    }
  }
}
