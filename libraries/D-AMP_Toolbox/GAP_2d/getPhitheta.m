function y = getPhitheta(theta,ydim)
%N = length(theta);
y_all = fwht(theta);
y = y_all(1:ydim,:);
end