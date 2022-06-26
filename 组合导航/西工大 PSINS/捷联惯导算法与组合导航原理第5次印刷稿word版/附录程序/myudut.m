function [U, D] = myudut(P)  % UD分解，P=U*diag(D)*U', U为上三角阵
    n = length(P);
    U = eye(n); D = zeros(n,1);  trPn = trace(P)/n*1e-40;
    for j=n:-1:1
        k = (j+1):n;
        D(j) = P(j,j) - (U(j,k).^2)*D(k);
        if D(j)<=trPn, continue; end
        for i=(j-1):-1:1
            U(i,j) = (P(i,j)-(U(i,k).*U(j,k))*D(k)) / D(j);
        end
end
