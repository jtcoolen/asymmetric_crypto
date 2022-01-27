#include <pollard.h>
#include <quadratic_sieve.h>
#include <stdio.h>

int main(void) {
  mpz_t fac, fac2, n;
  mpz_inits(fac, fac2, n, NULL);

  mpz_init_set_str(n, "52590354472497239257283147", 10);
  factor_rho_pollard(n, fac);
  gmp_printf("factor_rho_pollard(%Zu)=%Zu\n", n, fac);
  mpz_init_set_str(n, "52590354464570687296135717939981", 10);
  factor_rho_pollard(n, fac);
  gmp_printf("factor_rho_pollard(%Zu)=%Zu\n", n, fac);

  mpz_init_set_str(n, "1042387", 10);
  unsigned long long P = 50;
  unsigned long long A = 500;

  quadratic_sieve(n, P, A, fac, fac2);
  gmp_printf("quadratic_sieve(%Zu)=(%Zu, %Zu)\n", n, fac, fac2);

  mpz_clears(fac, fac2, n, NULL);

  return 0;
}
