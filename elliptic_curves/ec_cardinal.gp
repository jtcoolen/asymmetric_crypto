ellisoncurve2(E, P) = (P[2]^2 + E[1] * P[1] * P[2] + E[3] * P[2]) == (P[1]^3 + E[2] * P[1]^2 + E[4] * P[1] + E[5]);

\\ fonctionne sur un corps premier ou fini
ellcard_naive(E, q) = {
  my(c = 1, i, j, x, y, g); \\ point à l'infini
  qprime = isprime(q);
  if(qprime,
    g = znprimroot(q);
  );
  if(!qprime, g = ffprimroot(ffgen(q)));
  \\ x!=0 et y!=0
  x = g;
  for(i = 1, q - 1,
    y = g;
    for(j = 1, q - 1,
      if(ellisoncurve2(E, [x, y]), c++);
      y *= g;
    );
    x *= g;
  );
  \\ x==0 et y==0
  if(E[5] == 0, c++); \\ [0,0] est sur la courbe
  x = g;
  \\y==0 et x!=0
  for(i = 1, q - 1,
    if(0 == x^3 + E[4]*x + E[5], c++);
    x = x * g;
  ); 
  \\x==0 et y!=0
  y = g;
  for(i = 1, q - 1,
    if(y^2 == E[5], c++);
      y = y * g;
  );
  c
}

\\ fonctionne sur une courbe elliptique définie sur un corps premier
\\ ne fonctionne qu'en caractéristique différente de 2 et 3 !
ellcard_legendre(E) = {
  my(c = 1, i); \\ point à l'infini
  for(i = 0, E.p - 1, c += kronecker(lift(lift(x^3 + E.a4 * x + E.a6)), E.p) + 1); c
}

ellcard2(E) = {
  if(E.p == 2 || E.p == 3, return(ellcard_naive(E)));
  ellcard_legendre(E)
}

\\ test (corps premier)
p = randomprime(2^10);
a = Mod(2, p);

E  = ellinit([a^4, a^6], a);
print("E  = ", E );
print("ellcard(E ) = ", ellcard(E ));
print("ellcard_naive(E ) = ", ellcard_naive(E , p));
print("ellcard_legendre(E ) = ", ellcard_legendre(E ));

\\ test (corps fini)
n = 5^4;
g = ffprimroot(ffgen(n));
E = ellinit([g^3, g^4]);
print("E = ", E);
print("ellcard(E) = ", ellcard(E));
print("ellcard_naive(E) = ", ellcard_naive(E, n));