function dpos = deltapos(pos)
    cl = cos(pos(:,1)); Re = 6378137;
    dpos = [[(pos(:,1)-pos(1,1)),(pos(:,2)-pos(1,2)).*cl]*Re,pos(:,3)-pos(1,3)];