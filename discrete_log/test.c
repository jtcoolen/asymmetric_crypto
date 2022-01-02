#include <discrete_log.h>
#include <stdio.h>

int main(void) {
  mpz_t order, gen, h, log;
  mpz_inits(order, gen, h, log, NULL);

  mpz_set_ui(order, 113);
  mpz_set_ui(gen, 3);
  mpz_set_ui(h, 57);

  babystep_giantstep(order, gen, h, log);

  printf("log=");
  mpz_out_str(stdout, 10, log);

  mpz_set_ui(h, 229);     // 229);
  mpz_set_ui(gen, 2);     // 2);
  mpz_set_ui(order, 383); // 383);
  printf("\nPollard's Rho (O(1) memory)\n");
  pollard_rho(order, gen, h, log);
  printf("log=");
  mpz_out_str(stdout, 10, log);
  return 0;
}