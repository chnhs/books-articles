function qnb = qdelphi(qpb, phi)  % 四元数真实值=四元数计算值-失准角误差
    qnb = qmul(rv2q(phi), qpb);