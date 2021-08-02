function H = measurement_entropy(y, total_pixels)

y = (y(:)).^2;

nbin = round(mean(y));
frequency = hist(y,nbin);

p = frequency/sum(frequency);
p(p==0) = 1;
H = -sum(p.*log2(p));

H = H*numel(y)/total_pixels;

  
