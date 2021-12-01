#include <finite_field.h>
#include <stdlib.h>
#include <stdio.h>

int main(void) {
  uint64_t b[4] = {1,2,3,4};
  uint64_t b1[3] = {1, 2};
  struct prime_field_element *p = prime_field_element_from_array(b, 4, 3);
  struct prime_field_element *q = prime_field_element_from_array(b1, 2, 3);
  prime_field_element_print(p);
  printf("\n");
  prime_field_element_print(q);
  struct prime_field_element *pq = prime_field_element_multiplication(p, q);
  printf("\n");
  prime_field_element_print(pq);
  printf("\n");
  struct prime_field_element *rem, *quot;
  prime_field_element_division(p, q, &rem, &quot);
  prime_field_element_print(rem);
  printf("\n");
  prime_field_element_print(quot);
  printf("\n");
  uint64_t b2[14] = {0, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1};
  struct prime_field_element *r = prime_field_element_from_array(b2, 14, 3);
  uint64_t b3[13] = {0, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1};
  struct prime_field_element *s = prime_field_element_from_array(b3, 13, 3);
  struct prime_field_element *gcd_pr = prime_field_element_gcd(p, r);
  //prime_field_element_print(gcd_pr);
  printf("\n");
  struct finite_field_element * pr = malloc(sizeof(struct finite_field_element));
  pr->poly = prime_field_element_copy(p);
  pr->mod = prime_field_element_copy(r);

  struct finite_field_element * sr = malloc(sizeof(struct finite_field_element));
  sr->poly = prime_field_element_copy(s);
  sr->mod = prime_field_element_copy(r);

  struct finite_field_element *res = finite_field_element_multiplication(pr, sr);
  prime_field_element_print(res->poly);
  printf("\n");
  struct prime_field_element *u, *v;
  struct prime_field_element *gcd = prime_field_element_gcd_extended(p, r, &u, &v);
  printf("\n");
  prime_field_element_print(gcd);
  printf("\n");
  prime_field_element_print(u);
  printf("\n");
  prime_field_element_print(v);

  return 0;
}