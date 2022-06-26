function test_L_7
    n = 3; m = 2; l = 2;
    [Phi, Gamma, H, Q, R, P0] = rndmodel(n, m, l);
    % (1) 标准KF
    P10 = Phi*P0*Phi'+Gamma*Q*Gamma';
    P1 = P10 - P10*H'*(H*P10*H'+R)^-1*H*P10;
    % (2) 平方根KF
    sQ = mychol(Q); sR = mychol(R); Delta0 = mychol(P0);
    Delta1 = SRKF(Delta0, Phi, Gamma, sQ, H, sR);  
    errSRKF = P1 - Delta1*Delta1',
    % (3) 平方根信息KF
    iPhi = Phi^-1;  isQ = sQ^-1; isR = sR^-1; S0 = mychol(P0^-1);
    S1 = SRIKF(S0, iPhi, Gamma, isQ, H, isR);
    errSRIKF = P1 - (S1*S1')^-1,
    % (4) 奇异值分解KF
    [U, Lambda] = svd(P0);  Lambda = sqrt(Lambda);
    [U1, Lambda1] = SVDKF(U, Lambda, Phi, Gamma, sQ, H, isR);
    errSVD = P1 - U1*Lambda1^2*U1',
    % (5) UD分解KF
    [U, D] = myudut(P0);
    [U, D] = UDKF(U, D, Phi, Gamma, diag(Q), H(1,:), R(1,1), 'TM');
    for k=2:length(R)
        [U, D] = UDKF(U, D, Phi, Gamma, diag(Q), H(k,:), R(k,k), 'M');
    end
    errUD = P1 - U*diag(D)*U',
    
function [Phi, Gamma, H, Q, R, P0] = rndmodel(n, m, l)  % 随机系统模型
    Phi = randn(n);  Gamma = randn(n,l);  H = randn(m,n);
    Q = diag(randn(l,1))^2;  R = diag(randn(m,1))^2;
    P0 = randn(n); P0 = P0'*P0;
    
function Delta1 = SRKF(Delta0, Phi, Gamma, sQ, H, sR)  % 平方根滤波
    [q, Delta] = myqr([Phi*Delta0, Gamma*sQ]');  Delta = Delta';
    [q, rho] = myqr([H*Delta, sR]'); rho = rho';
    Delta1 = Delta*(eye(length(Delta0))-Delta'*H'*(rho*rho'+sR*rho')^-1*H*Delta);
        
function S = SRIKF(S0, iPhi, Gamma, isQ, H, isR)  % 平方根信息滤波
    [q, rho] = myqr([S0'*iPhi*Gamma; isQ]); rho = rho';
    S = iPhi'*S0*(eye(length(S0))-...
                  S0'*iPhi*Gamma*(rho*rho'+isQ*rho')^-1*Gamma'*iPhi'*S0);
    [q, S] = myqr([S'; isR'*H]); S = S';
    
function [U, Lambda] = SVDKF(U, Lambda, Phi, Gamma, sQ, H, isR)  % 奇异值分解滤波
    [U, Lambda] = svd([Phi*U*Lambda, Gamma*sQ]);
    [U, Lambda] = svd([U*diag(1./diag(Lambda)), H'*isR]);
    Lambda = diag(1./diag(Lambda));

function [U, D] = UDKF(U, D, Phi, Gamma, Q, H, R, TM)  % UD滤波
    n = length(U);
    if ~isempty(strfind(upper(TM),'T'))
        W = [Phi*U, Gamma];   D1 = [D; Q];  % D,Q为向量
        for j=n:-1:1   % 时间更新
            D(j) = (W(j,:).*W(j,:))*D1;
            for i=(j-1):-1:1
                U(i,j) = (W(i,:).*W(j,:))*D1/D(j);
                W(i,:) = W(i,:) - U(i,j)*W(j,:);
            end
        end
    end
    if ~isempty(strfind(upper(TM),'M'))
        f = (H*U)';  g = D.*f;  afa = f'*g+R;
        for j=n:-1:1   % 量测更新
            afa0 = afa - f(j)*g(j); lambda = -f(j)/afa0;
            D(j) = afa0/afa*D(j);   afa = afa0;
            for i=(j-1):-1:1
                s = (i+1):(j-1);
                U(i,j) = U(i,j) + lambda*(g(i)+U(i,s)*g(s));
            end
        end
    end
