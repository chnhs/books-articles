function shplot(n, m, cs)
% Spherical harmonics plot, with n-degree m-order.
    dtheta = pi/20;    theta = 0:dtheta/2:pi;    lambda = 0:dtheta:2*pi;
    the1 = repmat(theta, length(lambda), 1);
    lam1 = repmat(lambda', 1, length(theta));
    r = ones(size(the1));
    x = r.*sin(the1).*cos(lam1);  y = r.*sin(the1).*sin(lam1);  z = r.*cos(the1);
 
    P = legendre(n, cos(theta), 'norm');
    if cs==1  %  cs=1 for cos, cs=0 for sin.
        c = cos(m*lambda)'*P(m+1,:);
        str = sprintf('cos(%d\\it\\lambda\\rm)x\\itP\\rm_{%d}^{%d}(cos\\it\\theta)', m,n,m);
    else
        c = sin(m*lambda)'*P(m+1,:);
        str = sprintf('sin(%d\\it\\lambda\\rm)x\\itP\\rm_{%d}^{%d}(cos\\it\\theta\\rm)',...
                      m,n,m);
    end
    surf(x,y,z,c); axis equal; title(str);
