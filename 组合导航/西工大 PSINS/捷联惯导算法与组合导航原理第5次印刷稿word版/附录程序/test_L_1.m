gvar
x = linspace(-1,1);
for n=0:5
    P = legendre(n, x);
    for k=0:n
        subplot(2,3,k+1); xlabel('x'); ylabel(sprintf('P_n^%d(x)',k));
        if n==0, w=4; elseif n==1, w=1; else w=2; end
        hold on; grid on; plot(x, (-1)^k*P(k+1,:), lsc(n+1,:), 'linewidth',w);
        xlim([-1.03 1.03]);
    end
end
subplot(2,3,1); ylim([-1.1 1.1]);  legend('n=0','n=1','n=2','n=3','n=4','n=5');
