function qpb = qaddphi(qnb, phi)  % ��Ԫ������ֵ=��Ԫ����ʵֵ+ʧ׼�����
    qpb = qmul(rv2q(-phi),qnb);
