function [wm, vm] = av2imu(att, vn, pos, ts)  % 由姿态、速度和位置生成惯性传感器信息
    wm0 = zeros(3,1); vm0 = wm0; I33 = eye(3);
    wm = att(2:end,:); vm = wm;
    for k=2:length(att)
        eth = earth((pos(k-1,:)+pos(k,:))'/2, (vn(k-1,:)+vn(k,:))'/2);
        qbb = qmul(qmul(qconj(a2qua(att(k-1,:))),rv2q(eth.wnin*ts)),a2qua(att(k,:)));
        phim = q2rv(qbb);
        wm1 = (I33+askew(1/12*wm0))\phim;
        dvnsf = vn(k,:)'-vn(k-1,:)'-eth.gcc*ts;  Cnb0 = a2mat(att(k-1,:)');
        vm1 = (I33+1/2*askew(1/6*wm0+wm1))\...
              (Cnb0'*(I33+askew(eth.wnin*ts/2))*dvnsf-1/12*cross(vm0,wm1));
        wm(k-1,:) = wm1';  vm(k-1,:) = vm1;  wm0 = wm1; vm0 = vm1;
    end