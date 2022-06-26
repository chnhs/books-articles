clear all
n = 4; cs = 1; % cs=1 for cos, cs=0 for sin.
for m=0:n
    subplot(fix(n/6)+1, min(6,n+1), m+1);  shplot(n,m,cs);
end
