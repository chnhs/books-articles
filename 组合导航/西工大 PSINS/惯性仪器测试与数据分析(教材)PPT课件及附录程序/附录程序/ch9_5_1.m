function test_ekf
% 《惯性仪器测试与数据分析》EKF滤波仿真, Yan Gongmin, 2012-08-22
    Ts = 0.2; % 采样周期
    t = 10; % 仿真时间
    len = floor(t/Ts); % 仿真步数
    g = 9.8; rho = 0.05; % 重力, 阻尼系数
    % 真实轨迹模拟
    dax = 0.15; day = 0.1;  % 系统噪声
    dr = 10; dafa = 0.01; % 量测噪声
    X = zeros(len,4); X(1,:) = [0, 30, 200, 0]; % 状态模拟的初值
    Y(1,:) = hhh(X(1,:)')' + [dr, dafa].*randn(1,2);
    for k=2:len
        X(k,:) = fff(X(k-1,:)', rho, g, Ts)' + [0, dax*randn(1), 0, day*randn(1)];
        Y(k,:) = hhh(X(k,:)')' + [dr, dafa].*randn(1,2);
    end
    % EKF 滤波
    Qk = diag([0; dax; 0; day])^2; Rk = diag([dr; dafa])^2;
    Xk = X(1,:)'; Pk = 100*eye(4);
    for k=1:len
        [fX, Fk] = fff(Xk, rho, g, Ts);   
        [hX, Hk] = hhh(fX);
        [Xk, Pk, Kk] = ekf(Fk, Qk, Pk, Hk, Rk, Y(k,:)', fX, hX);
        X_est(k,:) = Xk';
    end
    plot(X(:,1),X(:,3),'-b', Y(:,1).*sin(Y(:,2)),Y(:,1).*cos(Y(:,2)),'*', X_est(:,1),X_est(:,3),'+r')
    xlabel('\itx'); ylabel('\ity');
    legend('real', 'measurement', 'EKF estimated'); grid on

function [fX, JF] = fff(X, rho, g, Ts) % 系统状态非线性函数
    x = X(1); vx = X(2); y = X(3); vy = X(4); 
    fX = X+[vx; -rho*vx^2; vy; rho*vy^2-g]*Ts;
    JF = diag([1, 1-2*rho*vx*Ts, 1, 1+2*rho*vy*Ts]); JF(1,2) = Ts; JF(3,4) = Ts;
    
function [hX, JH] = hhh(X) % 量测非线性函数
    x = X(1); y = X(3);
    r2 = x^2+y^2; r = sqrt(r2);
    hX = [r; atan(x/y)];
    JH = [x/r, 0, y/r, 0;   y/r2, 0, -x/r2, 0];

function [Xk, Pk, Kk] = ekf(Phi, Qk, Pk_1, Hk, Rk, Yk, fX, hX) % EKF 滤波函数
    Pkk_1 = Phi*Pk_1*Phi' + Qk;     Pxy = Pkk_1*Hk';    Pyy = Hk*Pxy + Rk;     
    Kk = Pxy*Pyy^-1;
    Xk = fX + Kk*(Yk-hX);
    Pk = Pkk_1 - Kk*Pyy*Kk';
