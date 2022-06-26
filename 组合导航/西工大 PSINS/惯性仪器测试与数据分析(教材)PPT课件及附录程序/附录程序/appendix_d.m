% 《惯性仪器测试与数据分析》惯导误差近似解析解的验证, by Yan Gongmin 2012-08-22
clear all
R = 6378137; g = 9.8; ug = 1e-6*g;
arcdeg = pi/180; arcmin = arcdeg/60; hur =3600; dph = arcdeg/hur; wie = 15.041067*dph;
L = 30*arcdeg;
Ts = 10;  T = 24*hur; t = [Ts:Ts:T]';
en = [1; 2; 3]*0.01*dph; Dn = [50; 100]*ug;
fi0 = [10;10;0]*arcmin; dv0 = [0; 0]; dLti0 = 0*arcmin; dLgi0 = 0*arcmin;
% fi0 = [-Dn(2)/g; Dn(1)/g; en(1)/(wie*cos(L))]; % 采用自对准
X0 = [fi0; dv0; dLti0; dLgi0]; U = [en; Dn; 0; 0];
% 解析解
fE0=X0(1); fN0=X0(2); fU0=X0(3); dvE0=X0(4); dvN0=X0(5); dLti0=X0(6); dLgi0=X0(7);
eE=U(1); eN=U(2); eU=U(3); DE=U(4); DN=U(5);
sL = sin(L); cL = cos(L); tL = tan(L); eL = sec(L); sL2 = sL^2; cL2 = cL^2;
wN = wie*cL; wU = wie*sL; wf = wie*sL; 
ws = sqrt(g/R); ws = sqrt(ws^2+wf^2); V1 = R*ws;
ss = sin(ws*t); cs = cos(ws*t); se = sin(wie*t); ce = cos(wie*t); sf = sin(wf*t); cf = cos(wf*t);
c11=-DE/g*cs.*sf; c12=-DN/g*(1-cs.*cf); c13=-dLti0*wie/ws*sL*ss.*sf;  % fE
    c14=fE0*cs.*cf; c15=fN0*cs.*sf; c16=-fU0*wie/ws*cL*ss.*cf;
    c17=eE/ws*ss.*cf; c18=eN/ws*ss.*sf; c19=0;
    c110=dvE0/V1*ss.*sf; c111=-dvN0/V1*ss.*cf;
c21=DE/g*(1-cs.*cf); c22=-DN/g*cs.*sf; c23=-dLti0*wie/ws*sL*ss.*cf;  % fN
    c24=-fE0*cs.*sf; c25=fN0*cs.*cf; c26=fU0*wie/ws*cL*ss.*sf;
    c27=-eE/ws*ss.*sf; c28=eN/ws*ss.*cf; c29=0; 
    c210=dvE0/V1*ss.*cf; c211=dvN0/V1*ss.*sf;
c31=DE/g*tL*(1-cs.*cf); c32=-DN/g*tL*cs.*sf; c33=dLti0*eL*(se-wie/ws*sL2*ss.*cf);  % fU
    c34=fE0*eL*(se-sL*cs.*sf); c35=fN0*tL*(cs.*cf-ce); c36=fU0*(ce+wie/ws*sL*ss.*sf);
    c37=eE*eL*((1-ce)/wie-sL/ws*ss.*sf); c38=-eN*tL*(se/wie-ss.*cf/ws); c39=eU/wie*se;
    c310=dvE0/V1*tL*ss.*cf; c311=dvN0/V1*tL*ss.*sf;
c41=DE/g*V1*ss.*cf; c42=DN/g*V1*ss.*sf; c43=dLti0*R*wU*(ce-cs.*cf);  % dvE
    c44=fE0*V1*ss.*sf; c45=-fN0*V1*ss.*cf; c46=fU0*R*wN*(cs.*sf-sL*se);
    c47=eE*R*(sL*se-cs.*sf); c48=eN*R*(cs.*cf-cL2-sL2*ce); c49=-eU*R*cL*(sL*(1-ce)-wie/ws*ss.*sf);
    c410=dvE0*cs.*cf; c411=dvN0*cs.*sf;
c51=-DE/g*V1*ss.*sf; c52=DN/g*V1*ss.*cf; c53=dLti0*R*wie*(sL*cs.*sf-se);  % dvN
    c54=fE0*V1*ss.*cf; c55=fN0*V1*ss.*sf; c56=fU0*R*wN*(cs.*cf-ce);
    c57=eE*R*(ce-cs.*cf); c58=eN*R*(sL*se-cs.*sf); c59=-eU*R*cL*(se-wie/ws*ss.*cf);
    c510=-dvE0*cs.*sf; c511=dvN0*cs.*cf;
c61=DE/g*cs.*sf; c62=DN/g*(1-cs.*cf); c63=dLti0*(ce+wie/ws*sL*ss.*sf);  % dLti
    c64=fE0*(ce-cs.*cf); c65=fN0*(sL*se-cs.*sf); c66=-fU0*cL*(se-wie/ws*ss.*cf);
    c67=eE*(se/wie-ss.*cf/ws); c68=eN*(sL/wie*(1-ce)-ss.*sf/ws); c69=-eU/wie*cL*(1-ce);
    c610=-dvE0/V1*ss.*sf; c611=dvN0/V1*ss.*cf;
c71=DE/g*eL*(1-cs.*cf); c72=-DN/g*eL*cs.*sf; c73=dLti0*tL*(se-wie/ws*ss.*cf);  % dLgi
    c74=fE0*eL*(sL*se-cs.*sf); c75=fN0*eL*(cs.*cf-cL2-sL2*ce); c76=-fU0*(sL*(1-ce)-wie/ws*ss.*sf);
    c77=eE*eL*(sL/wie*(1-ce)-ss.*sf/ws); c78=-eN*(cL*t+sL*tL/wie*se-eL/ws*ss.*cf); c79=-eU*sL*(t-se/wie);
    c710=dvE0*eL/V1*ss.*cf; c711=dvN0*eL/V1*ss.*sf; c712=dLgi0;
fE   = c11+c12+c13 +c14+c15+c16 +c17+c18+c19 +c110+c111;
fN   = c21+c22+c23 +c24+c25+c26 +c27+c28+c29 +c210+c211;
fU   = c31+c32+c33 +c34+c35+c36 +c37+c38+c39 +c310+c311;
dvE  = c41+c42+c43 +c44+c45+c46 +c47+c48+c49 +c410+c411;
dvN  = c51+c52+c53 +c54+c55+c56 +c57+c58+c59 +c510+c511;
dLti = c61+c62+c63 +c64+c65+c66 +c67+c68+c69 +c610+c611;
dLgi = c71+c72+c73 +c74+c75+c76 +c77+c78+c79 +c710+c711+c712;
X = [fE, fN, fU, dvE, dvN, dLti, dLgi];
% 数值解
F = [ 0   wU -wN  0    -1/R   0   0
     -wU  0   0   1/R   0    -wU  0  
      wN  0   0   tL/R  0     wN  0  
      0  -g   0   0     2*wU  0   0  
      g   0   0  -2*wU  0     0   0  
      0   0   0   0     1/R   0   0  
      0   0   0   eL/R  0     0   0 ];
[Fk, Bk] = c2d(F, eye(size(F)), Ts); Uk = Bk*U;  % 离散化
Xk = X0; XX = X;
for k=1:length(t)
    Xk = Fk*Xk+Uk;
    XX(k,:) = Xk';
end
figure(1)
subplot(411), plot(t/3600,[X(:,1:2),XX(:,1:2)]/arcmin); ylabel('\it\phi_E ,\phi_N / \prime'); grid
subplot(412), plot(t/3600,[X(:,3),XX(:,3)]/arcmin); ylabel('\it\phi_U / \prime'); grid
subplot(413), plot(t/3600,[X(:,4:5),XX(:,4:5)]); ylabel('\it\deltav_E ,\deltav_N / m/s'); grid
subplot(414), plot(t/3600,[X(:,6:7),XX(:,6:7)]/arcmin); ylabel('\it\deltaL ,\delta\lambda / \prime'); grid
xlabel('\itt \rm/ hur');

    return;
    figure(2)
    subplot(411), plot(t/3600,[X(:,1:2)-XX(:,1:2)]/arcmin); ylabel('arcmin'); grid
    subplot(412), plot(t/3600,[X(:,3)-XX(:,3)]/arcmin); ylabel('arcmin'); grid
    subplot(413), plot(t/3600,[X(:,4:5)-XX(:,4:5)]); ylabel('m/s'); grid
    subplot(414), plot(t/3600,[X(:,6:7)-XX(:,6:7)]/arcmin); ylabel('arcmin'); grid
    xlabel('t / hur')

    % 验证表D-1 8-10行
    c81=0; c82=0; c83=dLti0*ce;
        c84=fE0*ce; c85=fN0*sL.*se; c86=-fU0*cL.*se;
        c87=eE/wie*se; c88=eN/wie*sL*(1-ce); c89=-eU/wie*cL*(1-ce);
    c91=0; c92=0; c93=-dLti0*sL.*se;
        c94=-fE0*sL*se; c95=fN0*(cL2+sL2*ce); c96=fU0*sL*cL*(1-ce);
        c97=-eE/wie*sL*(1-ce); c98=eN*(cL2*t+sL2/wie*se); c99=eU*sL*cL*(t-se/wie); c912=-dLgi0*cL;
    c101=0; c102=0; c103=dLti0*cL.*se;
        c104=fE0*cL*se; c105=fN0*sL*cL*(1-ce); c106=fU0*(sL2+cL2*ce);
        c107=eE/wie*cL*(1-ce); c108=eN*sL*cL*(t-se/wie); c109=eU*(sL2*t+cL2/wie*se); c1012=-dLgi0*sL;
    fEdLti = c81+c82+c83 +c84+c85+c86 +c87+c88+c89;
    fNdLgi = c91+c92+c93 +c94+c95+c96 +c97+c98+c99 +c912;
    fUdLgi = c101+c102+c103 +c104+c105+c106 +c107+c108+c109 +c1012;
    apr = [ dLti0+fE0+fN0*sL*wie*t-fU0*cL*wie*t+eE*t, ... 
            -dLti0*sL*wie*t-fE0*sL*wie*t+fN0+eN*t-dLgi0*cL, ...
            dLti0*cL*wie*t+fE0*cL*wie*t+fU0+eU*t-dLgi0*sL ];
    figure(3)
    hold off, plot(t/3600,[[fE+dLti,fN-dLgi*cL,fU-dLgi*sL],[fEdLti,fNdLgi,fUdLgi]]/arcmin); ylabel('arcmin'); grid
    hold on, plot(t/3600,apr/arcmin); ylabel('arcmin'); grid
    % 验证表D-2
    apr0 = [ c11+c12+c14+c15, c21+c22+c24+c25, c31+c32+c34+c35, ...
             c41+c42+c44+c45, c51+c52+c54+c55, c61+c62+c64+c65, c71+c72+c74+c75 ];
    apr1 = [ c16+c17, c26+c27, c36+c37, c46+c47, c56+c57, c66+c67, c76+c77 ];
    II = ones(size(t));
    apr2 = [ -DN/g*II, DE/g*II, DE/g*tL*(1-ce)-DN/g*eL*se, 0*II, 0*II, ...
             DE/g*sL*se+DN/g*(1-ce), DE/g*sL*tL*(1-ce)-DN/g*tL*se ];
    figure(4)
    subplot(411), plot(t/3600,[XX(:,1:2),apr2(:,1:2)]/arcmin); ylabel('arcmin'); grid
    subplot(412), plot(t/3600,[XX(:,3),apr2(:,3)]/arcmin); ylabel('arcmin'); grid
    subplot(413), plot(t/3600,[XX(:,4:5),apr2(:,4:5)]); ylabel('m/s'); grid
    subplot(414), plot(t/3600,[XX(:,6:7),apr2(:,6:7)]/arcmin); ylabel('arcmin'); grid
    