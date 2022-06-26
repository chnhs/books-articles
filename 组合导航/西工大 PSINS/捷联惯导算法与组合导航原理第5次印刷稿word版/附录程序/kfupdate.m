function kf = kfupdate(kf, Zk, TimeMeasBoth)
    if nargin==1,         TimeMeasBoth = 'T';
    elseif nargin==2,     TimeMeasBoth = 'B';    end
    if TimeMeasBoth=='T' || TimeMeasBoth=='B'     % 时间更新
        kf.Xkk_1 = kf.Phikk_1*kf.Xk;
        kf.Pkk_1 = kf.Phikk_1*kf.Pk*kf.Phikk_1' + kf.Gammak*kf.Qk*kf.Gammak';
    else % TimeMeasBoth=='M'
        kf.Xkk_1 = kf.Xk;
        kf.Pkk_1 = kf.Pk; 
    end
    if TimeMeasBoth=='M' || TimeMeasBoth=='B'     % 量测更新
        kf.PXZkk_1 = kf.Pkk_1*kf.Hk';
        kf.PZkk_1 = kf.Hk*kf.PXZkk_1 + kf.Rk;
        kf.Kk = kf.PXZkk_1/kf.PZkk_1;
        kf.Xk = kf.Xkk_1 + kf.Kk*(Zk-kf.Hk*kf.Xkk_1);
        kf.Pk = kf.Pkk_1 - kf.Kk*kf.PZkk_1*kf.Kk';
    else % TimeMeasBoth=='T'
        kf.Xk = kf.Xkk_1;
        kf.Pk = kf.Pkk_1; 
    end
	kf.Pk = (kf.Pk+kf.Pk')/2;   % P阵对称化
