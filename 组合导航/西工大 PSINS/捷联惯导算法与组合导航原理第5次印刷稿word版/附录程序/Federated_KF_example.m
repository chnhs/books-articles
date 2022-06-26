function Federated_KF_example
% 问题描述参见《捷联惯导算法与组合导航原理》练习题48. by YanGongmin
    Ts = 0.1;  N = 1000;
    q = 0.01;  r = [10; 5];
    Ft = [0 1; 0 0];                     
    Hk = [1 0; 1 0];
    Phikk_1 = eye(2) + Ft*Ts;
    Gammak_1 = [0; 1];
    Qk = q*Ts;
    Rk = diag(r)^2;
    X0 = [1000; -10];
    Xk = X0;
    Xkk = zeros(N, 2); Zkk = zeros(N, 2);
    for k = 1:N
        Xk = Phikk_1*Xk + Gammak_1*sqrt(Qk)*randn(1);
        Zk = Hk*Xk + r.*randn(2,1);
        Xkk(k,:) = Xk';   Zkk(k,:) = Zk';
    end
    figure
    subplot(211), plot(Xkk(:,2)); xlabel('k / 0.1s'); ylabel('velocity / m/s'); grid on; title('trajectory simulation')
    legend('falling velocity');
    subplot(212), plot([Xkk(:,1),Zkk]); xlabel('k / 0.1s'); ylabel('height / m'); grid on
    legend('real height', 'height meas1', 'height meas2');
    %% standard KF
    Xk0 = X0+[10;3].*randn(2,1);
    Pk0 = diag([100; 10]);
    Xk = Xk0; Pk = Pk0;
    KF_res = zeros(N, 4);
    for k = 1:N
        [Xk, Pk] = Kalman(Phikk_1, Gammak_1, Qk, Xk, Pk, Hk, Rk, Zkk(k,:)');
        KF_res(k,:) = [Xk; diag(Pk)]';
    end
    %% sequential KF
    Xk = Xk0; Pk = Pk0;
    SKF_res = zeros(N, 4);
    for k = 1:N
        [Xk, Pk] = Kalman(Phikk_1, Gammak_1, Qk, Xk, Pk, 1);
        [Xk, Pk] = Kalman(Xk, Pk, Hk(1,:), Rk(1,1), Zkk(k,1));
        [Xk, Pk] = Kalman(Xk, Pk, Hk(2,:), Rk(2,2), Zkk(k,2));
        SKF_res(k,:) = [Xk; diag(Pk)]';
    end
    %% federated KF
    beta = rand(3,1); beta = beta/sum(beta);
    beta1 = beta(1); beta2 = beta(2); betam = beta(3);
    Xk = Xk0; Pk = Pk0;
    FKF_res = zeros(N, 4);
    for k = 1:N
        Xk1 = Xk; Xk2 = Xk; Xkm = Xk;
        Pk1 = Pk; Pk2 = Pk; Pkm = Pk;
        [Xk1, Pk1] = Kalman(Phikk_1, Gammak_1, Qk, Xk1, Pk1, Hk(1,:), Rk(1,1), Zkk(k,1), beta1);
        [Xk2, Pk2] = Kalman(Phikk_1, Gammak_1, Qk, Xk2, Pk2, Hk(2,:), Rk(2,2), Zkk(k,2), beta2);
        [Xkm, Pkm] = Kalman(Phikk_1, Gammak_1, Qk, Xkm, Pkm, betam);
%         [Xkm, Pkm] = Kalman(Phikk_1, Gammak_1, Qk, Xkm, Pkm/betam, 1);
        Pk = ((Pk1)^-1 + (Pk2)^-1 + (Pkm)^-1)^-1 ;
        Xk = Pk * ((Pk1)^-1*Xk1 + (Pk2)^-1*Xk2 + (Pkm)^-1*Xkm);
        FKF_res(k,:) = [Xk; diag(Pk)];
    end
    figure
    subplot(211), plot([KF_res(:,2)-Xkk(:,2), SKF_res(:,2)-Xkk(:,2), FKF_res(:,2)-Xkk(:,2)]); title('estimated errors')
    xlabel('k / 0.1s'); ylabel('velocity err / m/s'); grid on; legend('standard KF', 'sequential KF', 'federated KF');
    subplot(212), plot([KF_res(:,1)-Xkk(:,1), SKF_res(:,1)-Xkk(:,1), FKF_res(:,1)-Xkk(:,1)]);
    xlabel('k / 0.1s'); ylabel('height err / m'); grid on; legend('standard KF', 'sequential KF', 'federated KF');
    figure
    subplot(211), plot(sqrt([KF_res(:,4), SKF_res(:,4), FKF_res(:,4)])); title('sqrt(Pk)')
    xlabel('k / 0.1s'); ylabel('velocity err / m/s'); grid on; legend('standard KF', 'sequential KF', 'federated KF');
    subplot(212), plot(sqrt([KF_res(:,3), SKF_res(:,3), FKF_res(:,3)]));
    xlabel('k / 0.1s'); ylabel('height err / m'); grid on; legend('standard KF', 'sequential KF', 'federated KF');

function [Xk, Pk] = Kalman(Phikk_1, Gammak_1, Qk_1, Xk_1, Pk_1, Hk, Rk, Zk, beta)
    if nargin==5  % measure update: [Xk, Pk] = Kalman(Xkk_1, Pkk_1, Hk, Rk, Zk)
        Xkk_1 = Phikk_1; Pkk_1 = Gammak_1; Hk = Qk_1; Rk = Xk_1; Zk = Pk_1;
        Kk = Pkk_1*Hk'*(Hk*Pkk_1*Hk' + Rk)^-1;
        Xk = Xkk_1 + Kk*(Zk - Hk*Xkk_1);
        Pk = (eye(length(Xkk_1)) - Kk*Hk)*Pkk_1;
    elseif nargin==6  % time update: [Xkk_1, Pkk_1] = Kalman(Phikk_1, Gammak_1, Qk_1, Pk_1, Xk_1, beta)
        beta = Hk;
        Xkk_1 = Phikk_1*Xk_1;
        Pkk_1 = Phikk_1*Pk_1*Phikk_1' + Gammak_1*Qk_1*Gammak_1';
        Pk = Pkk_1/beta; Xk = Xkk_1;
    else% time/measure update: [Xk, Pk] = Kalman(Phikk_1, Gammak_1, Qk_1, Xk_1, Pk_1, Hk, Rk, Zk, beta)
        Xkk_1 = Phikk_1*Xk_1;
        Pkk_1 = Phikk_1*Pk_1*Phikk_1' + Gammak_1*Qk_1*Gammak_1';
        if nargin==9,  Pkk_1 = Pkk_1/beta;  end
        Kk = Pkk_1*Hk'*(Hk*Pkk_1*Hk' + Rk)^-1;
        Xk = Xkk_1 + Kk*(Zk - Hk*Xkk_1);
        Pk = (eye(length(Xkk_1)) - Kk*Hk)*Pkk_1; 
    end
