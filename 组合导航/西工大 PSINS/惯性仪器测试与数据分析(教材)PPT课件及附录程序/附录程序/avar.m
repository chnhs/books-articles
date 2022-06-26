function [sigma, tau, Err] = avar(y0, tau0)
% ����Allan����
% ���룺y -- ����(һ�л�������), tau0 -- ��������
% �����sigma -- Allan����(���ٵ�λ������y����һ��), tau -- ȡ��ʱ��, Err -- �ٷֱȹ������
% ����: Yan Gong-min, 2012-08-22
% example: 
%     y = randn(100000,1) + 0.00001*[1:100000]';
%     [sigma, tau, Err] = avar(y, 0.1);
    N = length(y0);
    y = y0; NL = N;
    for k = 1:inf
        sigma(k,1) = sqrt(1/(2*(NL-1))*sum([y(2:NL)-y(1:NL-1)].^2));
        tau(k,1) = 2^(k-1)*tau0;
        Err(k,1) = 1/sqrt(2*(NL-1));
        NL = floor(NL/2);
        if NL<3
            break;
        end
        y = 1/2*(y(1:2:2*NL) + y(2:2:2*NL));  % ���鳤�ȼӱ�(���ݳ��ȼ���)
    end
    subplot(211), plot(tau0*[1:N], y0); grid
    xlabel('\itt \rm/ s'); ylabel('\ity');
    subplot(212), 
    loglog(tau, sigma, '-+', tau, [sigma.*(1+Err),sigma.*(1-Err)], 'r--'); grid
    xlabel('\itt \rm/ s'); ylabel('\it\sigma_A\rm( \tau )');
