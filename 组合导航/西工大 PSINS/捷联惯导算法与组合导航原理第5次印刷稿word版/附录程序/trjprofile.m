function [att, vn, pos] = trjprofile(att0, vn0, pos0, wat, ts)  % 生成轨迹姿态、速度和位置参数
    len = fix(sum(wat(:,5))/ts);
    att = zeros(len, 3); vn = att; pos = att;  kk=1;
    att(1,:) = att0'; vn(1,:) = vn0'; pos(1,:) = pos0';
    vb = a2mat(att0)'*vn0; vby = vb(2);   % 求纵向速度
    b = fir1(20, 0.01, 'low');  b = b/sum(b); x = repmat([att0;vby]',length(b),1); % 低通滤波器
    for m=1:size(wat,1);
        watk = wat(m,:);
        for tk=ts:ts:(watk(5)+ts/10)
            att0 = att0 + watk(1:3)'*ts;   vby = vby + watk(4)*ts;
            x = [x(2:end,:); [att0;vby]']; y = b*x;  % 低通滤波
            att(kk+1,:) = y(1:3);
            vn(kk+1,:) = (a2mat(att(kk+1,:)')*[0;y(4);0])';
            vn01 = (vn(kk,:)+vn(kk+1,:))/2;
            eth = earth(pos(kk,:)',vn01');
            pos(kk+1,:) = pos(kk,:) + [vn01(2)/eth.RMh;vn01(1)/eth.clRNh;vn01(3)]'*ts;
            kk = kk+1;
        end
    end
    att(kk:end,:) = []; vn(kk:end,:) = []; pos(kk:end,:) = [];
    