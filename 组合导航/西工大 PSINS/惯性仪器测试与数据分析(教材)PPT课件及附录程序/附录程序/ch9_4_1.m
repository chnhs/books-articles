function test_akf
% 《惯性仪器测试与数据分析》自适应Kalman滤波仿真, Yan Gongmin, 2012-08-22
    % 数据生成仿真，模型: xk = phi*xk_1 + wk, yk = h*xk + vk
    N = 1000;
    lq = 500; q0 = [-.5*ones(lq,1); .5*ones(N-lq,1)];
    lQ = 700; Q0 = [5*ones(lQ,1); 1*ones(N-lQ,1)].^2;
    lr = 600; r0 = [2*ones(lr,1); 0.2*ones(N-lr,1)];
    lR = 300; R0 = [5*ones(lR,1); 30*ones(N-lR,1)].^2;
    phi = .99; h = 1;
    x = 0; y = 0;
    for k=2:N
        x(k,1) = phi*x(k-1,1) + (q0(k-1)+sqrt(Q0(k-1))*randn(1));
        y(k,1) = h*x(k,1) + (r0(k)+sqrt(R0(k))*randn(1));
    end
    % Sage-Husa 自适应Kalman滤波
    xk = 0; Pk = 10; qk = 0; Qk = 0; rk = 0; Rk = 100; RA = 100;
    b = 0.98; beta = 1;
    for k=2:N
        beta = beta/(beta+b);
        [xk, Pk, qk, Qk, rk, Rk] = akf(xk, Pk, phi, h, y(k), qk, Qk, rk, Rk, beta);
%        qk = q0(k);  % 注释了则表示不自适应
        Qk = Q0(k); 
        rk = r0(k); 
%         Rk = R0(k); 
        RA = (1-beta)*RA + beta/2*(y(k)-y(k-1))^2; Rk = RA; % Allan方差求量测噪声
        res(k,:) = [xk, qk, Qk, rk, Rk];
    end
    yy = [x, y, res(:,1)];
    subplot(311), plot(1:5:N, yy(1:5:end,:));xlabel('\itk'); ylabel('\ity_k , x_k , x^\^_k'); grid
    subplot(323), plot([q0, res(:,2)]); xlabel('\itk'); ylabel('\itq_k'); grid
    subplot(324), plot([R0, res(:,5)]); xlabel('\itk'); ylabel('\itR_k'); grid
    subplot(325), plot([Q0, res(:,3)]); xlabel('\itk'); ylabel('\itQ_k'); grid
    subplot(326), plot([r0, res(:,4)]); xlabel('\itk'); ylabel('\itr_k'); grid

function [xk, Pk, qk, Qk, rk, Rk] = akf(xk_1, Pk_1, phi, h, yk, qk, Qk, rk, Rk, beta)
    Pxkk_1 = phi*Pk_1*phi'+Qk; Pxykk_1 = Pxkk_1*h'; Pykk_1 = h*Pxykk_1+Rk;
    Kk = Pxykk_1*Pykk_1^-1;
    xkk_1 = phi*xk_1+qk; ykk_1 = h*xkk_1+rk;
    eykk_1 = yk-ykk_1;  xk = xkk_1 + Kk*eykk_1;
    Pk = Pxkk_1 - Kk*Pykk_1*Kk';
    rk = (1-beta)*rk + beta*(yk-h*xkk_1);  % 噪声自适应
    Rk = (1-beta)*Rk + beta*(eykk_1*eykk_1'-h*Pxykk_1);
    qk = (1-beta)*qk + beta*(xk-phi*xk_1);
    Qk = (1-beta)*Qk + beta*(Kk*eykk_1*eykk_1'*Kk'+Pk-phi*Pk_1*phi');
