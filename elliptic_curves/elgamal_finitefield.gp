\\ Version d'ElGamal qui opère sur des courbes elliptiques définies sur des corps finis


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

decomp(a) = {
  my(s = 0);
  my(h = a);
  while(h % 2 == 0, h = (h >> 1); s++);
  [s, h]
}

\\ multiplication 2^r-aire modifiée d'un point
ellpointmul2rary(E, P, n, r) = {
  my(multiplesP = vector(2^r+1, i, [0, 0]));
  multiplesP[1] = P;
  multiplesP[2] = ellpointdup(E, P);

  \\ précalcul, j commence à 1 pour remplir toutes les cases d'indice impair
  for(j = 1, 2^(r - 1), multiplesP[2 * j + 1] = ellpointadd(E, multiplesP[2 * j - 1], multiplesP[2]));
  
  my(i = logint(n, 2^r));
  my(z, c, s, h);
  my(Q = oo);
  while(i != -1,
  
    c = (n >> (i * r)) % (2^r);
    if(c != 0, [s, h] = decomp(c));
    
    if(c != 0,
      Q = ellpointmulbin(E, Q, 2^(r - s));
      Q = ellpointadd(E, Q, multiplesP[h]));
      
    if(c == 0, s = r);
    
    Q = ellpointmulbin(E, Q, 2^s);
    i--;
  );
  Q
}

\\ Q: public key
\\ d: private key
ElGamal_gen_keys(E, P) = {
  my(d = random(ellcard(E) - 2) + 2);
  [ellpointmulbin(E, P, d), d]
}

ElGamal_encrypt(E, Q, P, msg) = {
  my(k = random(ellcard(E) - 2) + 2);
  [ellpointmulbin(E, P, k), ellpointadd(E, msg, ellpointmulbin(E, Q, k))]
}

ElGamal_decrypt(E, d, ciphertext) = {
  ellpointsub(E, ciphertext[2], ellpointmul2rary(E, ciphertext[1], d, 4))
}

test() = {
  for(i=0,5,
    n = randomprime(2^10)^2;
    g = ffprimroot(ffgen(n));
    Es = ellinit([g^3, g^4]);
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