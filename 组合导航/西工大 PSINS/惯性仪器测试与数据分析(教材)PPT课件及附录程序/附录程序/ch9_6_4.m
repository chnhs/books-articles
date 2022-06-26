function test_ukf
% �������������������ݷ�����UKF �˲�����, by Yan Gongmin 2012-08-22
    % ����ģ��
    len = 100;
    q = 0.1; r = .3; Pw = q^2; Pv = r^2;
    x = q*randn(1); y = x+r*randn(1);
    for k=2:len
        x(k,1) = fx(x(k-1)) + q*randn(1);
        y(k,1) = hx(x(k)) + r*randn(1);
    end
    % UKF�˲�
    Xk = 0;    Pk = 10.0;
    for k=1:len
        [Xk,Pk] = ukf(Xk, Pk, Pw, Pv, y(k));
        Xkk(k,1) = Xk;
    end
    plot([y, x, Xkk]); grid on, legend('����', '״̬��ֵ', '״̬��ֵ')
    xlabel('\itk'); ylabel('\ity_k , x_k , x^\^_k'); 
    
function y = fx(x)    % ״̬����
    y = sin(x);
    
function y = hx(x)   % ���ⷽ��
    if x>0
        y = x;
    else
        y = 2*x;
    end
    
function [Xk, Pk] = ukf(Xk_1, Pk_1, Pw, Pv, Yk)
    [Xkk_1, Pxx] = ut(Xk_1, Pk_1, @fx, 1e-3,2,0); % ״̬UT�任
    Pxx = Pxx + Pw;
    [Ykk_1, Pyy, Pxy] = ut(Xkk_1, Pxx, @hx, 1e-3,2,0); % ����UT�任
    Pyy = Pyy + Pv;
    Kk = Pxy*Pyy^-1; % �˲�
    Xk = Xkk_1+Kk*(Yk-Ykk_1);
    Pk = Pxx-Kk*Pyy*Kk';
    
function [y, Pyy, Pxy, X, Y] = ut(x, Pxx, hfx, alpha, beta, kappa) % UT �任
    % ͬ��¼D.8
    n = length(x);
    lambda = alpha^2*(n+kappa) - n;
    gamma = sqrt(n+lambda);
    Wm = [lambda/gamma^2; repmat(1/(2*gamma^2),2*n,1)];  
    Wc = [Wm(1)+(1-alpha^2+beta); Wm(2:end)];
    sPxx = gamma*chol(Pxx)';    % Choleskey ���Ƿֽ�
    xn = repmat(x,1,n); 
    X = [x, xn+sPxx, xn-sPxx];
    Y(:,1) = feval(hfx, X(:,1)); m=length(Y); y = Wm(1)*Y(:,1); Y = repmat(Y,1,2*n+1);
    for k=2:1:2*n+1     % �����Ա任�������ֵ
        Y(:,k) = feval(hfx, X(:,k));
        y = y + Wm(k)*Y(:,k);
    end
    Pyy = zeros(n); Pxy = zeros(n,m);
    for k=1:1:2*n+1
        yerr = Y(:,k)-y;
        Pyy = Pyy + Wc(k)*yerr*yerr';  % ���������
        xerr = X(:,k)-x;
        Pxy = Pxy + Wc(k)*xerr*yerr';  % �������Э������
    end
