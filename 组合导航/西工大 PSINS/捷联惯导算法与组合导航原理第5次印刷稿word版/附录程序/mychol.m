function A = mychol(P)  % 乔莱斯基分解，P=A*A', A为上三角阵
    n = length(P);  A = zeros(n);
    for j=n:-1:1
        A(j,j) = sqrt(P(j,j)-A(j,j+1:n)*A(j,j+1:n)');
        for i=(j-1):-1:1
            A(i,j) = (P(i,j)-A(i,j+1:n)*A(j,j+1:n)')/A(j,j);
        end
    end
