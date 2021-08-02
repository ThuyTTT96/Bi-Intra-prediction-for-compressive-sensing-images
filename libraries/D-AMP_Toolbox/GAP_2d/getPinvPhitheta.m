function x1 = getPinvPhitheta(y,theta,ydim)
N = length(theta);
y_all = fwht(theta);
y2 = y_all(1:ydim,:);
y2 = y-y2;
y2_all = [y2; zeros(length(theta)-ydim,1)];
x1 = fwht(y2_all);
end