gvar;
nn = 2; ts = 0.1; nts = nn*ts; % 子样数和采样时间
att0 = [0; 0; 30]*arcdeg; qnb0 = a2qua(att0);
vn0 = [0;0;0]; pos0 = [34*arcdeg; 108*arcdeg; 100];
qnb = qnb0;  vn = vn0;  pos = pos0; % 姿态、速度和位置初始化
eth = earth(pos, vn);
wm = qmulv(qconj(qnb),eth.wnie)*ts;  vm = qmulv(qconj(qnb),-eth.gn)*ts;
wm = repmat(wm', nn, 1); vm = repmat(vm', nn, 1);  % 仿真静态IMU数据
phi = [0.1; 0.2; 3]*arcmin;  qnb = qaddphi(qnb, phi);  % 失准角
eb = [0.01;0.015;0.02]*dph; web = [0.001;0.001;0.001]*dpsh;   % 陀螺常值零偏，角度随机游走
db = [80;90;100]*ug; wdb = [1;1;1]*ugpsHz;  % 加速度计常值偏值，速度随机游走
Qk = diag([web; wdb; zeros(9,1)])^2*nts;
rk = [[0.1;0.1;0.1];[[10;10]/Re;10]];  Rk = diag(rk)^2;
P0 = diag([[0.1;0.1;10]*arcdeg; [1;1;1]; [[10;10]/Re;10];
           [0.1;0.1;0.1]*dph; [80;90;100]*ug])^2;
Hk = [zeros(6,3),eye(6),zeros(6)];
kf = kfinit(Qk, Rk, P0, zeros(15), Hk);  % kf滤波器初始化
len = fix(3600/ts);  % 仿真时长
err = zeros(len, 10);  xkpk = zeros(len, 2*kf.n+1); kk = 1;  t = 0; % 记录导航结果
for k=1:nn:len
    t = t + nts;
    [wm1, vm1] = imuadderr(wm, vm, eb, web, db, wdb, ts);
    [qnb, vn, pos, eth] = insupdate(qnb, vn, pos, wm1, vm1, ts);
    kf.Phikk_1 = eye(15) + kfft15(eth, q2mat(qnb), sum(vm1,1)'/nts)*nts;
    kf = kfupdate(kf);
    if mod(t,1)<nts
        gps = [vn0; pos0] + rk.*randn(6,1);  % GPS速度位置仿真
        kf = kfupdate(kf, [vn;pos]-gps, 'M');
    end
    qnb = qdelphi(qnb,kf.Xk(1:3));  kf.Xk(1:3) = 0;  % 反馈
    vn = vn - kf.Xk(4:6);  kf.Xk(4:6) = 0;
    pos = pos - kf.Xk(7:9);  kf.Xk(7:9) = 0;
    err(kk,:) = [qq2phi(qnb,qnb0); vn-vn0; pos-pos0; t]';
    xkpk(kk,:) = [kf.Xk; diag(kf.Pk); t]; kk = kk+1;
    if mod(t,500)<nts,  disp(fix(t));  end  % 显示进度
end
err(kk:end,:) = [];  xkpk(kk:end,:) = [];  tt = err(:,end);
% 状态真值与估计效果对比图
msplot(321, tt, err(:,1:2)/arcmin, '\it\phi\rm / ( \prime )');
legend('\it\phi\rm_E', '\it\phi\rm_N'), 
msplot(322, tt, err(:,3)/arcmin, '\it\phi\rm_U\rm / ( \prime )');
msplot(323, tt, err(:,4:6), '\delta\itv^n\rm / ( m.s^{-1} )');
legend('\delta\itv\rm_E', '\delta\itv\rm_N', '\delta\itv\rm_U')
msplot(324, tt, [err(:,7)*Re,err(:,8)*Re*cos(pos(1)),err(:,9)], '\delta\itp\rm / m');
legend('\delta\itL', '\delta\it\lambda', '\delta\ith')
msplot(325, tt, xkpk(:,10:12)/dph, '\it\epsilon\rm / ( (\circ).h^{-1} )');
legend('\it\epsilon_x','\it\epsilon_y','\it\epsilon_z')
msplot(326, tt, xkpk(:,13:15)/ug, '\it\nabla\rm / \mu\itg');
legend('\it\nabla_x','\it\nabla_y','\it\nabla_z')
% 均方差收敛图
spk = sqrt(xkpk(:,16:end-1));
msplot(321, tt, spk(:,1:2)/arcmin, '\it\phi\rm / ( \prime )');
legend('\it\phi\rm_E', '\it\phi\rm_N'), 
msplot(322, tt, spk(:,3)/arcmin, '\it\phi\rm_U\rm / ( \prime )');
msplot(323, tt, spk(:,4:6), '\delta\itv^n\rm / ( m.s^{-1} )');
legend('\delta\itv\rm_E', '\delta\itv\rm_N', '\delta\itv\rm_U')
msplot(324, tt, [[spk(:,7),spk(:,8)*pos(1)]*Re,spk(:,9)], '\delta\itp\rm / m');
legend('\delta\itL', '\delta\it\lambda', '\delta\ith')
msplot(325, tt, spk(:,10:12)/dph, '\it\epsilon\rm / ( (\circ).h^{-1} )');
legend('\it\epsilon_x','\it\epsilon_y','\it\epsilon_z')
msplot(326, tt, spk(:,13:15)/ug, '\it\nabla\rm / \mu\itg');
legend('\it\nabla_x','\it\nabla_y','\it\nabla_z')
