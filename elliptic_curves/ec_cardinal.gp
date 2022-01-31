ellisoncurve2(E, P) = (P[2]^2 + E[1] * P[1] * P[2] + E[3] * P[2]) == (P[1]^3 + E[2] * P[1]^2 + E[4] * P[1] + E[5]);

ellcard_naive(E) = {
  my(c = 1, i, j); \\ point à l'infini
  for(i = 0, E.p - 1, for(j = 0, E.p - 1, if(ellisoncurve2(E,[i,j]), c++))); c
}

\\ ne fonctionne qu'en caractéristique différente de 2 et 3 !
ellcard_legendre(E) = {
  my(c = 1, i); \\ point à l'infini
  for(i = 0, E.p - 1, c += kronecker(lift(i^3 + E.a4 * i + E.a6), E.p) + 1); c
}

ellcard2(E) = {
  if(E.p == 2 || E.p == 3, return(ellcard_naive(E)));
  ellcard_legendre(E)
}

\\ test
p = randomprime(2^10)
a = Mod(2, p)

Es = ellinit([a^4, a^6], a);
print("Es = ", Es);
print("ellcard(Es) = ", ellcard(Es));
print("ellcard_naive(Es) = ", ellcard_naive(Es));
print("ellcard_legendre(Es) = ", ellcard_legendre(Es));