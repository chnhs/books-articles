% �������������������ݷ����������׷���, Yan Gongmin, 2012-08-22
% �źŷ���
fs = 200; %����Ƶ��
N = 1024; %��������
t = [0:N-1]*1/fs;  
xn = 1 + 5*sin(2*pi*20*t) + 3*sin(2*pi*25*t) + 0.1*randn(1,N); % ֱ��+�����ź�+������
subplot(121), plot(t,xn); grid, xlabel('t / s'); ylabel('x(t)');
% ֱ�ӷ�������(����ͼ��)
Xk = fft(xn,N);
S1 = abs(Xk).^2/N/fs;  S1(2:end) = S1(2:end)*2; % ���㵥�߹�����
subplot(122), hold off, 
N2 = floor(N/2);
semilogy([0:N2-1]*fs/N, S1(1:N2)); grid, xlabel('f / Hz'); ylabel('PSD');
% ��ӷ�������(��غ�����)
Rx = xcorr(xn,'biased'); S2 = xn;
for k = 1:N
    S2(k) = Rx(N) + 2*Rx(N+1:end)*cos((k-1)*2*pi/N*[1:N-1])';
end
S2 = S2/fs; S2(2:end) = S2(2:end)*2;
hold on, semilogy([0:N2-1]*fs/N, S2(1:N2),'r');
% Matlab���Դ���psd()����������(Welch�Ľ�����ͼ��)
S3 = psd(xn,N)/fs;  S3(2:end) = S3(2:end)*2;
hold on, semilogy([0:N2-1]*fs/N, S3(1:N2),'m');
