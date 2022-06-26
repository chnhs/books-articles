gvar
load egm2190;  % download egm2190.mat: https://pan.baidu.com/s/1bo3QLtX
N = 150;  % order
lat = linspace(-89.5,89.5)*arcdeg;  lon = 109*arcdeg;  res = zeros(length(lat),3);
for k=1:length(lat)
    gL = somigliana(lat(k), 0);
    gn = egm(GM, Re, wie, ff, egmc(1:N,1:N), egms(1:N,1:N), lat(k), lon, 0);
    res(k,:) = [-gn([2,1])'/norm(gn), norm(gn)-gL];
end
msplot(121, lat/arcdeg, res(:,1:2)/arcsec, '\itL\rm / ( \circ )', 'Deflection / ( \prime\prime )'); 
legend('\it\eta', '\it\xi'); xlim([-90,90]);
msplot(122, lat/arcdeg, res(:,3)/ug, '\itL\rm / ( \circ )', '\delta \itg\rm / \mug'); 
xlim([-90,90]);
