function vo = qmulv(q, vi)  % ��Ԫ����ʸ��������άʸ������Ԫ������任
    qi = [0;vi];
    qo = qmul(qmul(q,qi),qconj(q));
    vo = qo(2:4,1);
    % vo = q2mat(q)*vi;
