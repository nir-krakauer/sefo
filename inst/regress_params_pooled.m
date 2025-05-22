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

#t-distribution parameters (mu, sigma, nu) for a linear regression fit at Xstars
#under the hierarchical model Ys(:, i)|betas(:,i),sigmas(i) ~ N(Xs(:, :, i)*betas(:,i), sigmas(i)^2*I with betas(:,i) ~ N(beta0, diag(S0)). The hyperparameters beta0, S0 are estimated from the data using approximate maximum likelihood. (Partial pooling of the coefficient estimates across series.)
#sizes: Xs: n*k*m, Ys: n*m, Xstars: 1*k*m
#mu, sigma, nu: m*1
#syntax: [mu, sigma, nu, beta0, S0] = regress_params_pooled (Xs, Ys, Xstars)

#Example:
# n = 10; k = 3; m = 50;
# Xs = randn (n, k, m); beta00 = randn(k, 1); s00 = abs(randn(k, 1)); beta = beta00 + randn(k, m).*s00; Ys = squeeze(sum(Xs .* reshape(beta, [1 k m]), 2)) + randn(n, m); Xstars = randn(1, k, m);
# [mu, sigma, nu, beta0, S0] = regress_params_pooled (Xs, Ys, Xstars);
## The obtained beta0 and S0 should in this case approximate the beta00 and s00.^2 used to generate the input data, with the approximation getting better as the number of instances m increases.

## Author: Nir Krakauer <nkrakauer@ccny.cuny.edu>


function [mu, sigma, nu, beta0, S0] = regress_params_pooled (Xs, Ys, Xstars)

[n k m] = size (Xs);

bs = bs_sd = nan (k, m);
s2s = nan(m, 1);
for i = 1:m
  [~, ~, ~, beta_hat, s2, W] = regress_params (Xs(:, :, i), Ys(:, i));
  bs(:, i) = beta_hat;
  bs_var(:, i) = s2 * diag(W);
  s2s(i) = s2;
endfor
beta0 = sum(bs ./ bs_var, 2) ./ sum(1 ./ bs_var, 2);
S0 = max(mean(meansq(bs - beta0, 2) - bs_var, 2), eps);
L = spdiags(1 ./ sqrt(S0), 0, k, k);
[pl, nl] = size(L);
mu = sigma = nu = nan (m, 1);
for i = 1:m
  [X_s y_s L_generalized_inv beta_null] = linear_regularization_standardize(Xs(:, :, i), Ys(:, i), speye(n)/s2s(i), L, beta0);
  [beta_hat_s, lambda2, ~, SSE, T, W_s] = do_inversion_standard_aicc (X_s, y_s, 1);
  beta_hat = L_generalized_inv*beta_hat_s + (beta_null + beta0);
  W = L_generalized_inv*W_s*L_generalized_inv';

  Xstar = Xstars(:, :, i);
  mu(i) = Xstar*beta_hat;

  nu(i) = T - (nl - pl); #include unregularized fitted parameters

  s2 = SSE / nu(i);
  ss = sum(Xstar.* (W*Xstar')', 2);
  sigma(i) = sqrt(s2 * (1 + ss));
endfor


%{
n = 10; k = 3; m = 50;
Xs = randn (n, k, m); beta00 = randn(k, 1); s00 = abs(randn(k, 1)); beta = beta00 + randn(k, m).*s00; Ys = squeeze(sum(Xs .* reshape(beta, [1 k m]), 2)) + randn(n, m); Xstars = randn(1, k, m);
[mu, sigma, nu, beta0, S0] = regress_params_pooled (Xs, Ys, Xstars);
%}