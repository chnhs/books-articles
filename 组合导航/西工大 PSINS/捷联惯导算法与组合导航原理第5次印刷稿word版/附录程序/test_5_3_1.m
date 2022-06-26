Qk = 1; Rk = 2 ; P0 = 3; Phikk_1 = 1; Hk =1; Gammak = 1;
kf = kfinit(Qk, Rk, P0, Phikk_1, Hk, Gammak);
T = 5; Pkk_1 = zeros(T,1); Pk = Pkk_1; Kk = Pk;
for k=1:T
    kf = kfupdate(kf, 0);
    Pkk_1(k) = kf.Pkk_1; Pk(k) = kf.Pk; Kk(k) = kf.Kk; 
end
figure, plot(1:T, Pkk_1, '-.v', 0:T, [P0;Pk], '--^', 1:T, Kk, '-o', ...
             [0;reshape([1:T;1:T],2*T,1)], [P0;reshape([Pkk_1';Pk'],2*T,1)], ':'); 
xlabel('\itk'); ylabel('\itK_k , P_{k/k-1} , P_k'); ylim([0, 5]);
legend('\itP_{k/k-1}', '\itP_k', '\itK_k');

