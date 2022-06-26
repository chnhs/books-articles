gvar;
nn = 2; ts = 0.1; nts = nn*ts; % �������Ͳ���ʱ��
att = [1; 1; 30]*arcdeg; vn = [0;0;0]; pos = [34*arcdeg; 108*arcdeg; 100];
qnb = a2qua(att);  % ��̬���ٶȺ�λ�ó�ʼ��
eth = earth(pos, vn);
wm = qmulv(qconj(qnb),eth.wnie)*ts;  vm = qmulv(qconj(qnb),-eth.gn)*ts;
wm = repmat(wm', nn, 1); vm = repmat(vm', nn, 1);  % ���澲̬IMU����
phi = [0.1; 0.2; 3]*arcmin;  qnb = qaddphi(qnb, phi);
len = fix(3600/ts);  % ����ʱ��
avp = zeros(len, 10);  kk = 1;  t = 0; % ��¼������� [att, vn, pos, t]
for k=1:nn:len
    t = t + nts;
    [qnb, vn, pos] = insupdate(qnb, vn, pos, wm, vm, ts);  %vn(3) = 0;  pos(3) = 100;
    avp(kk,:) = [q2att(qnb); vn; pos; t]';  kk = kk+1;
    if mod(t,500)<nts,  disp(fix(t));  end  % ��ʾ����
end
avp(kk:end,:) = [];  tt = avp(:,end);
msplot(221, tt, avp(:,1:2)/arcdeg, 'Att / ( \circ )'); legend('\it\theta','\it\gamma')
msplot(222, tt, avp(:,3)/arcdeg, '\psi / ( \circ )');
msplot(223, tt, avp(:,4:6), 'Vel / ( m.s^{-1} )'); legend('\itv\rm_E', '\itv\rm_N', '\itv\rm_U')
msplot(224, tt, deltapos(avp(:,7:9)), '\DeltaPos / m');
       legend('\Delta\itL', '\Delta\it\lambda', '\Delta\ith')