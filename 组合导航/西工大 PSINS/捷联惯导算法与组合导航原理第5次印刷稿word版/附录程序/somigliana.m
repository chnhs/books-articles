function [gL, gLh, gN] = somigliana(lat, hgt)
global GM Re ff wie ge gp
    Rp = (1-ff)*Re;
    m = Re*wie^2/ge; beta = 5/2*m-ff; beta1 = (5*m*ff-ff^2)/8; beta2 = 2*GM/sqrt(Re*Rp)^3;
	sL = sin(lat); cL = cos(lat); sL2 = sL^2; cL2 = cL^2;
    gL = (Re*ge*cL2+Rp*gp*sL2)/sqrt(Re^2*cL2+Rp^2*sL2);  % Eq.(3.2.4)
    gLh = ge*(1+beta*sL2-beta1*sL2^2)-beta2*hgt; % Eq.(3.2.23)
    gN = 8.08e-9*hgt*sin(2*lat); % Eq.(3.2.24)