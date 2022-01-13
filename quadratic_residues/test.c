#include <quadratic_residues.h>
#include <stdio.h>

int main(void) {
  printf("square root of 302 mod 2081 = %ld", shanks_tonelli(302, 2081));
  return 0;
}
