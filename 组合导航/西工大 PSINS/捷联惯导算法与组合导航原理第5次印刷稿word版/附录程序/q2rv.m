function rv = q2rv(q) % �任��Ԫ��ת��Ϊ��Ч��תʸ��
	if q(1)<0,  q = -q;  end
    nmhalf = acos(q(1));  % ��Ч��תʸ��ģֵ��һ��
    if nmhalf>1e-20,   b = 2*nmhalf/sin(nmhalf);
    else             b = 2;    end
    rv = b*q(2:4);   % q = [ cos(|rv|/2); sin(|rv|/2)/|rv|*rv ];