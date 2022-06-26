function test_kf2
% 《惯性仪器测试与数据分析》Kalman滤波仿真, Yan Gongmin, 2012-08-22
    Phik = [0.95,0; 0 1];  Bk = eye(2);  Hk = [1.0 0];  % 系统结构参数
    q = diag([1; 1]); r = 3; Qk = q^2; Rk = r^2;  % 噪声参数
    len = 100; % 仿真步数
    % 随机误差模拟
    w = [q(1,1)*randn(len,1), q(2,2)*randn(len,1)]; v = r*randn(len,1);
    xk = zeros(len,2); yk = zeros(len,1);
    xk(1) = r*randn(1,1);
    for k=2:len
        xk(k,:) = (Phik*xk(k-1,:)' + Bk*w(k,:)')';
        yk(k) = Hk*xk(k,:)' + v(k);
    end
    % Kalman滤波估计
    Xk = [0;0]; Pxk = diag([10,10])^2;
    res = zeros(len,6);
    for k=1:len
        [Xk, Pxk, Kk] = kalman(Phik, Bk, Qk, Xk, Pxk, Hk, Rk, yk(k));
        res(k,:) = [Xk;sqrt(diag(Pxk));Kk]';
    end
    % 作图
    subplot(121), hold off, plot(xk(:,1),'-'), hold on, plot(res(:,1),'r:'); plot(res(:,3),'m-'); grid
        xlabel('\itk'); ylabel('\itX_{k1} , X^\^_{k1} , \surd P_{xk1}');
    subplot(122), hold off, plot(xk(:,2),'-'), hold on, plot(res(:,2),'r:'); plot(res(:,4),'m-'); grid
        xlabel('\itk'); ylabel('\itX_{k2} , X^\^_{k2} , \surd P_{xk2}');

function [Xk, Pxk, Kk] = kalman(Phikk_1, Bk, Qk, Xk_1, Pxk_1, Hk, Rk, Yk)
    Xkk_1 = Phikk_1*Xk_1;
    Pxkk_1 = Phikk_1*Pxk_1*Phikk_1' + Bk*Qk*Bk';
    Pxykk_1 = Pxkk_1*Hk';
    Pykk_1 = Hk*Pxykk_1 + Rk;
    Kk = Pxykk_1*Pykk_1^-1;
    Xk = Xkk_1 + Kk*(Yk-Hk*Xkk_1);
    Pxk = Pxkk_1 - Kk*Pykk_1*Kk';
