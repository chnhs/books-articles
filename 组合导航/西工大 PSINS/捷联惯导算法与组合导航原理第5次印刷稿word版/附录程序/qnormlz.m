function qnb = qnormlz(qnb)  % 四元数归一化
    nm = qnb'*qnb;
    if nm<1e-6,  qnb = [1; 0; 0; 0];
    else         qnb = qnb/sqrt(nm);    end