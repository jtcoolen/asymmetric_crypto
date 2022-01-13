#include <pollard.h>
#include <quadratic_sieve.h>
#include <stdio.h>

int main(void) {
  mpz_t fac, n;
  mpz_inits(fac, n, NULL);

  /*mpz_init_set_str(n, "52590354472497239257283147", 10);
  factor_rho_pollard(n, fac);
  gmp_printf("fac=%Zu\n", fac);
  mpz_init_set_str(n, "52590354464570687296135717939981", 10);
  factor_rho_pollard(n, fac);
  gmp_printf("fac=%Zu\n", fac);*/

  
  mpz_init_set_str(n, "1042387", 10);
  unsigned long long P = 50;
  unsigned long long A = 500;

  printf("sqrt(302) in Z/2081 Z %ld", shanks_tonelli(302, 2081));
  //quadratic_sieve(n, P, A);
    
  mpz_clears(fac, n, NULL);
  return 0;
}
