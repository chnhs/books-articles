function phi = qq2phi(qpb, qnb)  % ʧ׼�����=��Ԫ������ֵ-��Ԫ����ֵ
    qerr = qmul(qnb, qconj(qpb));
    phi = q2rv(qerr);
