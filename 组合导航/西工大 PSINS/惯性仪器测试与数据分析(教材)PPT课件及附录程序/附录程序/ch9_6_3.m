function test_UT_demo
% �������������������ݷ�������ά UT �任��ͼ��ʾ, by Yan Gongmin 2012-08-22
    % ��������
    mu = [0; 0];  % ��ֵ
    P11 = 1; P22 = 2; r = 0.42;   % rΪ���ϵ�� |r|<=1
    P12 = r*sqrt(P11*P22); P21 = P12;
    P0 = [P11,P12; P21,P22]; % ������
    subplot(121), hold off,
    plot(mu(1),mu(2), 'mx'); axis equal, hold on, ezplot(ezstring(mu,P0)); % ������Բ
    % ���ؿ������ӱ任
    x = mvnrnd(mu, P0, 500);  % �������ؿ�������
    plot(x(:,1), x(:,2), '.g'); grid on
    [mui, Pi] = parstat(x);
    plot(mui(1),mui(2), 'm*'); ezplot(ezstring(mui, Pi)); % ���ӷ�����Բ
    y = fx(x')';
    subplot(122), hold off, plot(y(:,1), y(:,2), '.g'); axis equal, grid on
    [muo, Po] = parstat(y);
    hold on, plot(muo(1),muo(2), 'm*');  ezplot(ezstring(muo, Po)); % ���ӷ�����Բ
    % �ƹ�����UT�任
    [y1, Pyy1, Pxy1, X1, Y1] = ut(mu, P0, @fx, 1, 0, 0);
    subplot(121), plot(X1(1,:),X1(2,:), 'ro');
    subplot(122), plot(Y1(1,:),Y1(2,:), 'ro');
    plot(y1(1),y1(2), 'rx'); ezplot(ezstring(y1, Pyy1));
     % �Ľ�������UT�任
    [y2, Pyy2, Pxy2, X2, Y2] = ut(mu, P0, @fx, 1e-3, 2, 0);
    subplot(121), plot(X2(1,:),X2(2,:), 'b^');
    subplot(122), plot(Y2(1,:),Y2(2,:), 'b^');
    plot(y2(1),y2(2), 'b+'); ezplot(ezstring(y2, Pyy2));
               
function z = fx(x)  % �Զ�������Ա任
    z = [(x(1,:)-1).*(x(2,:)-0.2); -(x(1,:)-1).^2];

function [mu, P] = parstat(x) % ��ά��������(��ֵ,����)ͳ��
    mu = mean(x,1);  % ��ֵ
    x = [x(:,1)-mu(1), x(:,2)-mu(2)];
    P = [x(:,1)'*x(:,1),x(:,1)'*x(:,2);x(:,2)'*x(:,1),x(:,2)'*x(:,2)]/length(x);
    
function ss = ezstring(mu, P) % ����ezplot�ַ���
    pp = P^-1;  
    ss = sprintf('%.2f*(x-%.2f)^2 + %.2f*(x-%.2f)*(y-%.2f) + %.2f*(y-%.2f)^2 - 1', ...
        pp(1,1), mu(1), 2*pp(1,2), mu(1), mu(2), pp(2,2), mu(2));
    
function [y, Pyy, Pxy, X, Y] = ut(x, Pxx, hfx, alpha, beta, kappa) % UT �任
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
