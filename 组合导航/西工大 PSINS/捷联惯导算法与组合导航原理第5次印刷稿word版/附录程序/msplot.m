function msplot(mnp, x, y, xstr, ystr)
    if mod(mnp,10)==1, figure; end   % 如果是第一幅小图，则新建一个figure
    subplot(mnp);
    if size(y,2)==2,      plot(x, y(:,1), '-', x, y(:,2), '-.');
    elseif size(y,2)==3   plot(x, y(:,1), '-', x, y(:,2), '-.', x, y(:,3), '--');
    else                  plot(x, y);     end
    % grid on;
    if nargin==4, ystr = xstr; xstr = '\itt\rm / s'; end  % 如果只输入一个字符串，则默认xlabel为时间
    xlabel(xstr); ylabel(ystr);