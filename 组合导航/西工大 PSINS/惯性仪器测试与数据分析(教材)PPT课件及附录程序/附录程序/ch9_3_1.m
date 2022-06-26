function test_kf
% 《惯性仪器测试与数据分析》Kalman滤波仿真, Yan Gongmin, 2012-08-22
    Phik = 0.95;  Bk = 1.0;  Hk = 1.0;  % 系统结构参数
    q = 1; r = 3; Qk = q^2; Rk = r^2;  % 噪声参数
    len = 100; % 仿真步数
    % 随机误差模拟
    w = q*randn(len,1); v = r*randn(len,1);
    xk = zeros(len,1); yk = zeros(len,1);
    xk(1) = r*randn(1,1);
    for k=2:len
        xk(k) = Phik*xk(k-1) + Bk*w(k);
        yk(k) = Hk*xk(k) + v(k);
    end
    % Kalman滤波估计
    Xk = 0; Pxk = 100*Rk/(Hk^2*Phik^2);
    for k=1:len
        [Xk, Pxk, Kk] = kalman(Phik, Bk, Qk, Xk, Pxk, Hk, Rk, yk(k));
        res(k,:) = [Xk,Pxk,Kk];
    end
    % 稳态滤波
    ss = [Hk^2*Phik^2  Hk^2*Bk^2*Qk+Rk-Phik^2*Rk  -Bk^2*Qk*Rk];
    Px = ( -ss(2) + sqrt(ss(2)^2-4*ss(1)*ss(3)) ) / (2*ss(1));
    K = Hk*(Phik^2*Px+Bk^2*Qk)/(Hk^2*(Phik^2*Px+Bk^2*Qk)+Rk); 
    G = (1-K*Hk)*Phik;
    Xk_IIR = filter(K, [1 -G], yk);
    % 作图
    subplot(121), hold off, plot(sqrt(res(:,2)),'-'), hold on, plot(res(:,3),'r:'); grid
        xlabel('\itk'); ylabel('\it\surd P_x_k , K_k');
    subplot(122), hold off, plot(yk,'x'), 
        hold on, plot(xk,'m:'); plot(res(:,1),'k'); plot(Xk_IIR,'r-.'); grid
        xlabel('\itk'); ylabel('\ity_k, x_k, x^\^_k, x^\^_k_,_I_I_R');

function [Xk, Pxk, Kk] = kalman(Phikk_1, Bk, Qk, Xk_1, Pxk_1, Hk, Rk, Yk)
    Xkk_1 = Phikk_1*Xk_1;
    Pxkk_1 = Phikk_1*Pxk_1*Phikk_1' + Bk*Qk*Bk';
    Pxykk_1 = Pxkk_1*Hk';
    Pykk_1 = Hk*Pxykk_1 + Rk;
    Kk = Pxykk_1*Pykk_1^-1;
    Xk = Xkk_1 + Kk*(Yk-Hk*Xkk_1);
    Pxk = Pxkk_1 - Kk*Pykk_1*Kk';
