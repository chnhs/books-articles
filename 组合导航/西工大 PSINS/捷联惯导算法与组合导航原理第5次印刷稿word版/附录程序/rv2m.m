function m = rv2m(rv)   % ��Ч��תʸ��ת��Ϊ�任����
	nm2 = rv'*rv;  % ��תʸ����ģ��
    if nm2<1.e-8   % ���ģ����С�������̩��չ��ǰ���������Ǻ���
        a = 1-nm2*(1/6-nm2/120); b = 0.5-nm2*(1/24-nm2/720);  % a->1, b->0.5
    else
        nm = sqrt(nm2);
        a = sin(nm)/nm;  b = (1-cos(nm))/nm2;
    end
    VX = askew(rv);
    m = eye(3) + a*VX + b*VX^2;