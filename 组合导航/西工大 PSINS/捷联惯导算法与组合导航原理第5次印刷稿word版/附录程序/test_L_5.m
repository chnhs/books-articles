gvar
L = 30*arcdeg; wN = wie*cos(L); wU = wie*sin(L); tL = tan(L); eL = sec(L);
Ts = 10;  T = 24*hur; t = (0:Ts:T-Ts)'/3600; len = length(t);
Ut = [-[0.01;0.0;0.0]*dph; [00;00]*ug; 0; 0];
X = zeros(7,len); X(:,1) = [[0;0;0]*arcmin; [0;0]; [0;0]/Re];
Ft = [ 0    wU  -wN   0     -1/Re   0   0
      -wU   0    0    1/Re   0     -wU  0  
       wN   0    0    tL/Re  0      wN  0  
       0   -g0   0    0      2*wU   0   0  
       g0   0    0   -2*wU   0      0   0  
       0    0    0    0      1/Re   0   0  
       0    0    0    eL/Re  0      0   0 ];
[Fk, Bk] = c2d(Ft, eye(size(Ft)), Ts); Uk = Bk*Ut;  % 离散化
for k=2:length(t)
    X(:,k) = Fk*X(:,k-1) + Uk;
%     X([4:7],k) = 0; % 纯失准角
%     X([1,3,5,6,7],k) = 0; % 东向通道
%     X([2,3,4,6,7],k) = 0; % 北向通道
%     X([3],k) = 0; % 水平通道
%     X([2,4],k) = 0; % 北向和方位通道
end
msplot(221, t, X(1:2,:)'/arcmin, '\itt\rm / h', '\it\phi\rm / ( \prime )');
legend('\it\phi\rm_E', '\it\phi\rm_N')
msplot(222, t, X(3,:)'/arcmin, '\itt\rm / h', '\it\phi\rm_U / ( \prime )');
msplot(223, t, X(4:5,:)', '\itt\rm / h', '\delta\itv^n\rm / ( m.s^{-1} )');
legend('\delta\itv\rm_E', '\delta\itv\rm_N')
msplot(224, t, X(6:7,:)'/arcmin, '\itt\rm / h', '\delta\itp\rm / ( \prime )');
legend('\delta\itL', '\delta\it\lambda')
