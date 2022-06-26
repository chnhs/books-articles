function [phim, dvbm] = cnscl(wm, vm)  % Ô²×¶Îó²î/»®´¬Îó²î²¹³¥
    cs = [  [2,    0,    0,    0,    0    ]/3
            [9,    27,   0,    0,    0    ]/20
            [54,   92,   214,  0,    0    ]/105
            [250,  525,  650,  1375, 0    ]/504
            [2315, 4558, 7296, 7834, 15797]/4620  ]; % 2-6×ÓÑù²¹³¥ÏµÊý
    wmm = sum(wm,1);  vmm = sum(vm,1);  dphim = zeros(1,3); scullm = zeros(1,3);
    n = size(wm, 1);  % ×ÓÑùÊý
	if n>1
        csw = cs(n-1,1:n-1)*wm(1:n-1,:); csv = cs(n-1,1:n-1)*vm(1:n-1,:);
        dphim = cross(csw,wm(n,:));  % Ô²×¶²¹³¥Á¿
        scullm = cross(csw,vm(n,:))+cross(csv,wm(n,:));  % »®´¬²¹³¥Á¿
    end
	phim = (wmm+dphim)';
	dvbm = (vmm+0.5*cross(wmm,vmm)+scullm)';