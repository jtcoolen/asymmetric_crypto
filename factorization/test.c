#include <pollard.h>


int main(void) {
  mpz_t fac, n;
  mpz_inits(fac, n, NULL);

  mpz_init_set_str(n, "52590354472497239257283147", 10);
  factor_rho_pollard(n, fac);
  gmp_printf("fac=%Zu\n", fac);
  mpz_init_set_str(n, "52590354464570687296135717939981", 10);
  factor_rho_pollard(n, fac);
  gmp_printf("fac=%Zu\n", fac);

  mpz_clears(fac, n, NULL);
  return 0;
}
