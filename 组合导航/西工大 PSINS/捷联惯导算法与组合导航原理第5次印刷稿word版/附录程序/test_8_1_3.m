gvar;  % ����ȫ�ֱ���
ts = 0.01;
att0 = [0;0;90]*arcdeg; vn0 = [0;0;0]; pos0 = [[34;108]*arcdeg;100];
%       ����������  ���������  ��λ������  ������ٶ�  ����ʱ��
wat = [ 0,         0,          0,         0,         10       %��ֹ
        0,         0,          0,         1,         10       %����
        0,         0,          0,         0,         10       %����
        5,         0,          0,         0,         4        %̧ͷ
        0,         0,          0,         0,         10       %����
       -5,         0,          0,         0,         4        %��ͷ
        0,         0,          0,         0,         10       %����
        0,         10,         0,         0,         1        %���
        0,         0,          9,         0,         10       %ת��
        0,        -10,         0,         0,         1        %���
        0,         0,          0,         0,         10       %����
        0,         0,          0,        -1,         10       %����
        0,         0,          0,         0,         10  ];   %��ֹ
wat(:,1:3) = wat(:,1:3)*arcdeg/1;  % ->deg/s
[att, vn, pos] = trjprofile(att0, vn0, pos0, wat, ts);
[wm, vm] = av2imu(att, vn, pos, ts);
tt = (0:length(att)-1)'*ts;
% �켣��ͼ
msplot(221, tt, att/arcdeg, 'Att / ( \circ )'); legend('\it\theta', '\it\gamma', '\it\psi')
msplot(222, tt, vn, 'Vel / ( m.s^{-1} )'); legend('\itv\rm_E', '\itv\rm_N', '\itv\rm_U')
msplot(223, tt, deltapos(pos), '\DeltaPos / m');
       legend('\Delta\itL', '\Delta\it\lambda', '\Delta\ith')
msplot(224, pos(:,2)/arcdeg, pos(:,1)/arcdeg, '\itL\rm / ( \circ )', '\it\lambda\rm / ( \circ )');
       hold on, plot(pos(1,2)/arcdeg, pos(1,1)/arcdeg, 'ro');
% ����������Ϣ��ͼ
msplot(121, tt(2:end), wm/ts/arcdeg, '\it\omega^b_{ib}\rm / ( (\circ).s^{-1} )');
       legend('\it\omega^b_{ibx}', '\it\omega^b_{iby}', '\it\omega^b_{ibz}');
msplot(122, tt(2:end), vm/ts, '\itf^b\rm_{sf} / ( m.s^{-2} )');
       legend('\itf^b\rm_{sf\itx}', '\itf^b\rm_{sf\ity}', '\itf^b\rm_{sf\itz}');

