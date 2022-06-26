gvar;  % 加载全局变量
ts = 0.01;
att0 = [0;0;90]*arcdeg; vn0 = [0;0;0]; pos0 = [[34;108]*arcdeg;100];
%       俯仰角速率  横滚角速率  方位角速率  纵向加速度  持续时间
wat = [ 0,         0,          0,         0,         10       %静止
        0,         0,          0,         1,         10       %加速
        0,         0,          0,         0,         10       %匀速
        5,         0,          0,         0,         4        %抬头
        0,         0,          0,         0,         10       %匀速
       -5,         0,          0,         0,         4        %低头
        0,         0,          0,         0,         10       %匀速
        0,         10,         0,         0,         1        %横滚
        0,         0,          9,         0,         10       %转弯
        0,        -10,         0,         0,         1        %横滚
        0,         0,          0,         0,         10       %匀速
        0,         0,          0,        -1,         10       %减速
        0,         0,          0,         0,         10  ];   %静止
wat(:,1:3) = wat(:,1:3)*arcdeg/1;  % ->deg/s
[att, vn, pos] = trjprofile(att0, vn0, pos0, wat, ts);
[wm, vm] = av2imu(att, vn, pos, ts);
tt = (0:length(att)-1)'*ts;
% 轨迹作图
msplot(221, tt, att/arcdeg, 'Att / ( \circ )'); legend('\it\theta', '\it\gamma', '\it\psi')
msplot(222, tt, vn, 'Vel / ( m.s^{-1} )'); legend('\itv\rm_E', '\itv\rm_N', '\itv\rm_U')
msplot(223, tt, deltapos(pos), '\DeltaPos / m');
       legend('\Delta\itL', '\Delta\it\lambda', '\Delta\ith')
msplot(224, pos(:,2)/arcdeg, pos(:,1)/arcdeg, '\itL\rm / ( \circ )', '\it\lambda\rm / ( \circ )');
       hold on, plot(pos(1,2)/arcdeg, pos(1,1)/arcdeg, 'ro');
% 惯性器件信息作图
msplot(121, tt(2:end), wm/ts/arcdeg, '\it\omega^b_{ib}\rm / ( (\circ).s^{-1} )');
       legend('\it\omega^b_{ibx}', '\it\omega^b_{iby}', '\it\omega^b_{ibz}');
msplot(122, tt(2:end), vm/ts, '\itf^b\rm_{sf} / ( m.s^{-2} )');
       legend('\itf^b\rm_{sf\itx}', '\itf^b\rm_{sf\ity}', '\itf^b\rm_{sf\itz}');

