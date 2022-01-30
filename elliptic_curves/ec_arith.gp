\\ Le point à l'infini est représenté par oo dans nos calculs
\\ On travaille en coordonnées affines
\\ Si l'on emploie des coordonnées projectives, alors un représentant du point à l'infini est (0,1,0)

ellpointneg(E, P) = {
  if(P == oo, return(P));
  [P[1], -E.a1 * P[1] - E.a3 - P[2]]
}

\\ addition en caractéristique différente de 2 et 3
ellpointadd(E, P1, P2) = {
  if(P1 == oo, return(P2));
  if(P2 == oo, return(P1));
  if(E.p == 2 || E.p == 3,
    s = elladd(E, P1, P2);
    if(s == [0], s = oo);
    return(s));
  P1 = Mod(P1,	Es.p);
  P2 = Mod(P2, Es.p);
  if(P1 == ellpointneg(E, P2), return(oo));
  if(P1 == P2,
    r = (3 * P1[1]^2 + E.a4) / (2 * P1[2]);
    x = r^2 - 2 * P1[1];
    return ([x, r * (P1[1] - x) - P1[2]]);
  );
  r = ((P2[2] - P1[2]) / (P2[1] - P1[1]));
  x = r^2 - P1[1] - P2[1];
  return ([x, r * (P1[1] - x) - P1[2]]);
}

ellpointsub(E, P1, P2) = ellpointadd(E, P1, ellpointneg(E, P2))

\\ duplication en caractéristique quelconque
ellpointdup(E, P) = {
  if(P == oo, return(P));
  if(P[2] == 0, return(ellpointadd(E, P, P))); \\ car il y a un risque que le dénominateur du gradient ci-dessous ne soit pas inversible
  \\ Gradient de la tangente au point P:
  tangentgradient = (3 * P[1]^2 + 2 * E.a2 * P[1] - E.a1 * P[2] + E.a4) / (2 * P[2] + E.a1 * P[1] + E.a3);
  x = tangentgradient^2 + tangentgradient * E.a1 - E.a2 - 2 * P[1];
  [x, -E.a1 * x - E.a3 - tangentgradient * x + tangentgradient * P[1] - P[2]]
}

\\ multiplication binaire d'un point
ellpointmulbin(E, P, n) = {
  my(Q = oo);
  my(k = logint(n, 2));
  while(k != -1, \\ pour visiter le cas k = 0
    Q = ellpointdup(E, Q);
    if((n >> k) % 2 == 1, Q = ellpointadd(E, P, Q));
    k--;
  );
  Q
}

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
  if(E.p == 2 || E.p == 3, return -1); \\ ne fonctionne qu'en caractéristique différente de 2 et 3
  my(c = 0);
  for(x = 0, E.p - 1,
    c += kronecker(lift(x^3 + E.a4 * x + E.a6), E.p) + 1
  );
  c++; \\ point à l'infini
  c
}

ellcard2(E) = {
  if(E.p == 2 || E.p == 3, return(ellcard_naive(E)));
  ellcard_legendre(E)
}

\\ Q: public key
\\ d: private key
ElGamal_gen_keys(E, P) = {
  my(d = random(ellcard(E) - 2) + 2);
  [ellpointmulbin(E, P, d), d]
}

ElGamal_encrypt(E, Q, P, msg) = {
  my(C = [oo, oo]);
  my(k=0);
  while(C[1] == oo || C[2] == oo,
    k = random(ellcard2(E) - 2) + 2;
    C = [ellpointmulbin(E, P, k), ellpointadd(E, msg, ellpointmulbin(E, Q, k))];
  );
  C
}

ElGamal_decrypt(E, d, ciphertext) = {
  ellpointsub(E, ciphertext[2], ellpointmulbin(E, ciphertext[1], d))
}

test() = {
  success = 1;
  while(success,
    p = randomprime(2^10);
    a = Mod(2, p);
    Es = ellinit([a^4, a^6], a);
    if(Es != [],
      print("Es = ", Es);
      P = [0]; \\ point à l'infini selon PARI/GP
      while(P == [0], P = random(Es));
      print("P = ", P);
      
      [Q, d] = ElGamal_gen_keys(Es, P);
      msg = [0];
      while(msg == [0], msg = random(Es));
      print("plaintext = ", msg);
      ciphertext = ElGamal_encrypt(Es, Q, P, msg);
      print("ciphertext = ", ciphertext);
      print("decrypted = ", ElGamal_decrypt(Es, d, ciphertext));
      success = (msg == ElGamal_decrypt(Es, d, ciphertext));
      print("Success: ", success)));
}

test()