% �������������������ݷ���������Ư�Ʒ���, Yan Gongmin, 2012-08-22
arcdeg = pi/180; hur = 3600; dph = arcdeg/hur; Hz = 1; % ���õ��ĵ�λ
eb = 0.1*dph*randn(1,1); % ��ֵƯ��
tauG = 50; beta = 1/tauG; R0 = 0.01*dph^2; %һ������ɷ���̵����ʱ���뷽��
q = 2*beta*R0; % �������������ϵ��������ʽ��9.1-16��
N2 = 0.0001*dph^2/Hz; fB = 400*Hz; % �۲������Ĺ����������
fs = 10; Ts = 1/fs; % ����Ƶ��,����
t = 600;  len = floor(t/Ts); % ����ʱ�䳤��
er = zeros(len,1); er(1) = sqrt(R0)*randn(1,1);
Phi = 1-beta*Ts; sQkr = sqrt(q*Ts); % һ������ɷ������ɢ��
for k=2:len
    er(k) = Phi*er(k-1) + sQkr * randn(1,1);
end
sQkg = sqrt(N2*fB); % �۲�����������
wg = sQkg*randn(len,1); 
subplot(121), plot([1:len]*Ts, [eb+er+wg,eb+er]/dph); grid on % ����ͼ
xlabel('\itt \rm/ s'); ylabel('\it\epsilon \rm/ (\circ)/h');
p1 = psd((er+wg)/dph,1024);    p2 = psd(er/dph,1024);
subplot(122), semilogy([0:512]*fs/1024, [p1/fs,p2/fs]), grid on % ������
xlabel('\itf \rm/ Hz'); ylabel('\itS\epsilon \rm/ ((\circ)/h)^2/Hz');
