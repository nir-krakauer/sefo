## Copyright (C) 2015-2016 Nir Krakauer
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; If not, see <http://www.gnu.org/licenses/>.

#t-distribution parameters (mu, sigma, nu) for a linear regression fit at Xstar
#and/or the estimated regression coefficient values (beta_hat), data error variance (s2), and coefficient covariance matrix (s2*W)
#syntax: [mu, sigma, nu, beta_hat, s2] = regress_params (X, Y[, Xstar])
#
#Example:
#n = 10; k = 8; p = 3; m = 2;
#X = randn(n, k); b = randn(k, p); Y = X*b + randn(n, p); Xstar = randn(m, k);
#[~, ~, ~, beta_hat, s2, W] = regress_params (X, Y);
#[mu, sigma, nu, beta_hat, s2] = regress_params (X, Y, Xstar);

## Author: Nir Krakauer <nkrakauer@ccny.cuny.edu>




function [mu, sigma, nu, beta_hat, s2, W] = regress_params (X, Y, Xstar)

n = size(Y, 1);
k = size(X, 2);
nu = n - k;

if isargout (6)
  W = pinv (X'*X);
endif

beta_hat = X \ Y; #or W*X'*Y
r = Y - X * beta_hat;
SSE = sum (r .^ 2); #residual sum of squares
s2 = SSE / nu; #same as (1/nu) * Y'*(speye(n) - X*W*X')*Y, but uses much less memory when n is large

if isargout (1)
  mu = Xstar*beta_hat;
endif

if isargout (2)
  #sumsq(Xstar / X, 2) and sum(Xstar.* (W*Xstar')', 2) [or sum(Xstar.* (pinv(X'*X)*Xstar')', 2))] are same as diag(Xstar*W*Xstar'), but more efficient when k is large
  mem_size = 1E8; #maximum size of array it's reasonable to temporarily create
  m = size(Xstar, 1);
  if n*m <= mem_size
    ss = sumsq(Xstar / X, 2);
  else
    if exist ("W", "var")
      ss = sum(Xstar.* (W*Xstar')', 2);
    else
      ss = sum(Xstar.* (pinv(X'*X)*Xstar')', 2);
    endif
  endif

  sigma = sqrt((1 + ss) * s2);
endif


%{
Xs: n*k*m
Ys: n*m
Xstars: 1*k*m
bs: k*m

[n k m] = size (Xs);

bs = bs_sd = nan (k, m);
s2s = nan(m, 1);
for i = 1:m
  [~, ~, ~, beta_hat, s2, W] = regress_params (Xs(:, :, i), Y(:, i));
  bs(:, m) = beta_hat;
  bs_var(:, m) = s2 * diag(W);
  s2s(m) = s2s;
endfor
b0 = sum(bs ./ bs_var, 2) ./ sum(1 ./ bs_var, 2);
S0 = max(mean(sumsq(bs - b0, 2) - bs_var, 2), 0);
L = spdiags(1 ./ sqrt(S0), 0, k, k);
for i = 1:m
  [X_s y_s L_generalized_inv beta_null] = linear_regularization_standardize(Xs(:, :, i), Y(:, i), speye(n)/s2s(m), L, b0);
  [beta_hat_s, lambda2, ~, SSE, T, W_s] = do_inversion_standard_aicc (X_s, y_s, 1);
  beta_hat = L_generalized_inv*beta_hat_s + (beta_null + beta0);
  W = L_generalized_inv*W_s*L_generalized_inv';

  mu(i) = Xstars(:, :, i)*beta_hat;

  [p, n] = size(L);
  nu(i) = T - (n - p); #include unregularized fitted parameters

  s2 = SSE / nu;
  ss = sum(Xstar.* (W*Xstar')', 2);
  sigma(i) = sqrt(s2 * (1 + ss));
endfor



%}
