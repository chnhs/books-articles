function att = q2att(qnb)  % 四元数转换为姿态角，以姿态阵作为中间变量
	att = m2att(q2mat(qnb));

    