function kf = kfinit(Qk, Rk, P0, Phikk_1, Hk, Gammak)
    [kf.m, kf.n] = size(Hk);
    kf.Qk = Qk; kf.Rk = Rk; kf.Pk = P0; kf.Xk = zeros(kf.n,1);
    kf.Phikk_1 = Phikk_1; kf.Hk = Hk;
    if nargin<6,  kf.Gammak = eye(kf.n);
    else          kf.Gammak = Gammak;   end

