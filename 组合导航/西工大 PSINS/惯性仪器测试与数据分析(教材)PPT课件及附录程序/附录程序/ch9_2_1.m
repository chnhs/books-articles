% 《惯性仪器测试与数据分析》陀螺漂移仿真, Yan Gongmin, 2012-08-22
arcdeg = pi/180; hur = 3600; dph = arcdeg/hur; Hz = 1; % 需用到的单位
eb = 0.1*dph*randn(1,1); % 常值漂移
tauG = 50; beta = 1/tauG; R0 = 0.01*dph^2; %一阶马尔可夫过程的相关时间与方差
q = 2*beta*R0; % 角速率随机游走系数，根据式（9.1-16）
N2 = 0.0001*dph^2/Hz; fB = 400*Hz; % 观测噪声的功率谱与带宽
fs = 10; Ts = 1/fs; % 采样频率,周期
t = 600;  len = floor(t/Ts); % 仿真时间长度
er = zeros(len,1); er(1) = sqrt(R0)*randn(1,1);
Phi = 1-beta*Ts; sQkr = sqrt(q*Ts); % 一阶马尔可夫过程离散化
for k=2:len
    er(k) = Phi*er(k-1) + sQkr * randn(1,1);
end
sQkg = sqrt(N2*fB); % 观测噪声均方差
wg = sQkg*randn(len,1); 
subplot(121), plot([1:len]*Ts, [eb+er+wg,eb+er]/dph); grid on % 序列图
xlabel('\itt \rm/ s'); ylabel('\it\epsilon \rm/ (\circ)/h');
p1 = psd((er+wg)/dph,1024);    p2 = psd(er/dph,1024);
subplot(122), semilogy([0:512]*fs/1024, [p1/fs,p2/fs]), grid on % 功率谱
xlabel('\itf \rm/ Hz'); ylabel('\itS\epsilon \rm/ ((\circ)/h)^2/Hz');
