function phi = qq2phi(qpb, qnb)  % 失准角误差=四元数计算值-四元数真值
    qerr = qmul(qnb, qconj(qpb));
    phi = q2rv(qerr);
