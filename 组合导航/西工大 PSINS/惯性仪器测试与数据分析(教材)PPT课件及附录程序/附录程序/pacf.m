function phi = pacf(rho)
% ����ARMA���̵�ƫ�����ϵ������(Levinson-Durbin recursion)
% ���� rho    --- �����ϵ��������������������rho(0)=1��
% ��� phi_kk --- ƫ�����ϵ��������������������phi_00=1��
% ����: Yan Gong-min, 2012-08-22
% example:
%     N = 1000;  k = floor(sqrt(N));
%     z1 = 0.9; z2 = -0.8; a1 = z1+z2; a2 = -z1*z2;
%     xn = filter(1, [1;-a1;-a2], randn(N,1)); % AR(2)
%     rho = xcorr(xn,'coeff'); rho = rho(N+1:N+k);
%     phi = pacf(rho(1:k));
%     subplot(211), plot(xn), grid
%     subplot(212), hold off, parcorr(xn,k), hold on
%     plot([0:k],[1,1;rho(1:k),phi(1:k)]); % ����acf��������������
    N = length(rho);
    phi(1) = rho(1); phi_kj = rho(1); % ��ֵ
    for k=1:N-1
        phi_kj(k+1,1) = (rho(k+1)-rho(k:-1:1)'*phi_kj) / (1-rho(1:k)'*phi_kj);
        phi_kj(1:k) = phi_kj(1:k) - phi_kj(end)*phi_kj(end-1:-1:1);
        phi(k+1,1) = phi_kj(end); % ���豣�����һ��ֵ
    end
