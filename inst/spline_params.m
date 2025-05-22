#t-distribution parameters (mu, sigma, nu) for smoothing spline fit at xpred
#with automatic choice of smoothing parameter, or it may be specified
#syntax: [mu, sigma, nu] = spline_params (x, y, xpred, p = [])

function [mu, sigma, nu, p] = spline_params (x, y, xpred, p = [])

n = numel (x);

if isempty (p)
  [~, p] = csaps_sel (x, y); #spline smoothing parameter selection
endif

[mu, p, sigma2, unc_yi, df] = csaps (x, y, p, xpred);

sigma = sqrt(sigma2 + unc_yi^2);

nu = n - df;
