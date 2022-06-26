function Ft = kfft15(eth, Cnb, fb)
global g0 ff
	tl = eth.tl; secl = 1/eth.cl; L = eth.pos(1); h = eth.pos(3); 
    f_RMh = 1/eth.RMh; f_RNh = 1/eth.RNh; f_clRNh = 1/eth.clRNh;
    f_RMh2 = f_RMh*f_RMh;  f_RNh2 = f_RNh*f_RNh;
    vE_clRNh = eth.vn(1)*f_clRNh; vE_RNh2 = eth.vn(1)*f_RNh2; vN_RMh2 = eth.vn(2)*f_RMh2;
    Mp1 = [ 0,           0, 0;
           -eth.wnie(3), 0, 0;
            eth.wnie(2), 0, 0 ];
    Mp2 = [ 0,             0,  vN_RMh2;
            0,             0, -vE_RNh2;
            vE_clRNh*secl, 0, -vE_RNh2*tl];
    beta = 5.27094e-3; beta1 = (2*beta+ff)*ff/8; beta2 = 3.086e-6; beta3 = 8.08e-9;
    Mp3 = [0,0,0; -2*beta3*h,0,-beta3*sin(2*L); -g0*(beta-4*beta1*cos(2*L))*sin(2*L),0,beta2];
    Maa = askew(-eth.wnin);
    Mav = [ 0,       -f_RMh, 0;
            f_RNh,    0,     0;
            f_RNh*tl, 0,     0 ];
    Map = Mp1+Mp2;
    Mva = askew(Cnb*fb);
    Mvv = askew(eth.vn)*Mav - askew(eth.wnien);
    Mvp = askew(eth.vn)*(Mp1+Map)+Mp3;
    Mpv = [ 0,       f_RMh, 0;
            f_clRNh, 0,     0;
            0,       0,     1 ];
    Mpp = [ 0,           0, -vN_RMh2;
            vE_clRNh*tl, 0, -vE_RNh2*secl;
            0,           0,  0 ];
    O33 = zeros(3);
	%%    phi   dvn   dpos   eb     db
	Ft = [ Maa    Mav    Map    -Cnb     O33 
           Mva    Mvv    Mvp     O33     Cnb 
           O33    Mpv    Mpp     O33     O33
           zeros(6,15) ];