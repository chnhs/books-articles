function [Q, R] = myqr(A)  % QR�ֽ⣬A=Q*R, ����Q'*Q=I��RΪ��������
    [m, n] = size(A);
    if n>m,  error('n must not less than m.'); end
    R = zeros(n);
    for i=1:n
        R(i,i) = sqrt(A(:,i)'*A(:,i));
        A(:,i) = A(:,i)/R(i,i);
        j = i+1:n;
        R(i,j) = A(:,i)'*A(:,j);
        A(:,j) = A(:,j)-A(:,i)*R(i,j);
    end
    Q = A;
