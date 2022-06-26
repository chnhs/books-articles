global GM Re ff wie ge gp g0 ug arcdeg arcmin arcsec hur dph dpsh ugpsHz lsc % Global VARables
GM = 3.986004415e14;  Re = 6.378136998405e6;  wie = 7.2921151467e-5;  % WGS-84 model
ff = 1/298.257223563; ee = sqrt(2*ff-ff^2); e2 = ee^2; Rp = (1-ff)*Re;
ge = 9.780325333434361; gp = 9.832184935381024; g0 = ge; ug = g0*1e-6;   % gravity, ug
arcdeg = pi/180; arcmin = arcdeg/60; arcsec = arcmin/60;   % angle unit
hur = 3600; dph = arcdeg/hur; dpsh = arcdeg/sqrt(hur); % hour, deg/hour, deg/sqrt(hour)
ugpsHz = ug/sqrt(1); % ug/sqrt(Hz)
lsc = [' -k';' -b';' -r';'-.m';'--g';' :c']; % line shape & color
