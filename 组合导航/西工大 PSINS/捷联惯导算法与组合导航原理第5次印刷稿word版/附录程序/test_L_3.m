gvar; 
C = zeros(11);  S = C;
C(1:2:11,1) = [1;                       -0.108262982131*10^-2/sqrt(5); % WGS-84 coefficients
        0.237091120053*10^-5/sqrt(9);   -0.608346498882*10^-8/sqrt(13); 
        0.142681087920*10^-10/sqrt(17); -0.121439275882*10^-13/sqrt(21) ];
lat = linspace(0,89.99)*arcdeg; lon = 109*arcdeg; hgt = [0, 5000, 10000, 20000];
for k=1:length(hgt)
    res = zeros(length(lat),5);
    for k1=1:length(lat)
        [gL, gLh, gN] = somigliana(lat(k1), hgt(k));
        gn = egm(GM, Re, wie, ff, C, S, lat(k1), lon, hgt(k));
        res(k1,:) = [gn', gLh, gN];
    end
    subplot(121);  plot(lat/arcdeg, -res(:,2)./abs(res(:,3))/arcsec, 'linewidth',2), hold on
                   plot(lat/arcdeg, res(:,5)./abs(res(:,3))/arcsec, 'r', 'linewidth',1),
    subplot(122);  plot(lat/arcdeg, (-res(:,3)-res(:,4))/ug, lsc(k,:)), hold on
end
subplot(121); xlabel('\itL\rm / ( \circ )'); ylabel('\it\xi\rm / ( \prime\prime )');
subplot(122); xlabel('\itL\rm / ( \circ )'); ylabel('\delta \itg\rm / \mug');
legend('\ith\rm = 0 m', '\ith\rm = 5,000 m', '\ith\rm = 10,000 m', '\ith\rm = 20,000 m');

