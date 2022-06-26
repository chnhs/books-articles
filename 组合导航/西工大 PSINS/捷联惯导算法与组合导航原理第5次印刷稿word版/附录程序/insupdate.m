function [qnb, vn, pos, eth] = insupdate(qnb, vn, pos, wm, vm, ts)  % �����ߵ���ֵ�����㷨
    nn = size(wm,1);  nts = nn*ts;
    [phim, dvbm] = cnscl(wm, vm);  % Բ׶���/��������
    eth = earth(pos, vn);  % ������ز�������
    vn1 = vn + rv2m(-eth.wnin*nts/2)*qmulv(qnb,dvbm) + eth.gcc*nts;  % �ٶȸ���
    vn = (vn+vn1)/2;
    pos = pos + [vn(2)/eth.RMh;vn(1)/eth.clRNh;vn(3)]*nts;  vn = vn1;  % λ�ø���
    qnb = qmul(rv2q(-eth.wnin*nts), qmul(qnb, rv2q(phim)));  % ��̬����
    qnb = qnormlz(qnb);