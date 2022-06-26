function test_L_6
    P = randn(5); P = P'*P;
    A = mychol(P);
    [U, D] = myudut(P);
    [Q, R] = myqr(P);
    errChol = norm(A*A' - P, 1),
    errUdut = norm(U*diag(D)*U' - P, 1),
    errQr = norm(Q*R - P, 1),
    