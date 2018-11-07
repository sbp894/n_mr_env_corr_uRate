function y=demean_scale(x, A)
y= x-min(x);
y=y/max(y)*A;