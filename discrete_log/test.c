#include <discrete_log.h>
#include <stdio.h>

int main(void) {
  mpz_t order, gen, h, log;
  mpz_inits(order, gen, h, log, NULL);

  mpz_set_ui(order, 113);
  mpz_set_ui(gen, 3);
  mpz_set_ui(h, 57);

  babystep_giantstep(order, gen, h, log);

  printf("\nlog=");
  mpz_out_str(stdout, 10, log);

  return 0;
}