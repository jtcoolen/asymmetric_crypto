ellcard_naive(E) = {
  my(c = 0);
  for(i = 0, E.p - 1,
    for(j = 0, E.p - 1,
      if(ellisoncurve(E, [i, j]), c++)
    )
  );
  c++; \\ point à l'infini
  c
}

ellcard_legendre(E) = {
  my(c = 0);
  for(x = 0, E.p - 1,
    c += kronecker(lift(x^3 + E.a4 * x + E.a6), Es.p) + 1
  );
  c++; \\ point à l'infini
  c
}

p = randomprime(2^10)
a = Mod(2, p)

Es = ellinit([a^4, a^6], a);
print("Es = ", Es);
print("ellcard(Es) = ", ellcard(Es));
print("ellcard_naive(Es) = ", ellcard_naive(Es));
print("ellcard_legendre(Es) = ", ellcard_legendre(Es));