## Copyright (C) 2015 Nir Krakauer
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
#with regularization (shrinkage toward zero) of the fitted coefficients
#syntax: [mu, sigma, nu] = regress_reg_gen_aicc_params (X, y, Xstar, C, L, beta0)
#assumed distributions: y ~ N(X*beta, sigma^2*C), L*(beta-beta0) ~ N(0, sigma^2*lambda^2*I), ystar ~ N(Xstar*beta, sigma^2); posterior t distribution accounts for uncertainty in beta and sigma but uses the AICC point estimate for lambda
#[mu, sigma, nu] = regress_reg_gen_aicc_params (X, y, Xstar, C, L, beta0);
#generalization of regress_reg_aicc_params

## Author: Nir Krakauer <nkrakauer@ccny.cuny.edu>


function [mu, sigma, nu] = regress_reg_gen_aicc_params (X, y, Xstar, C, L, beta0)

[X_s y_s L_generalized_inv beta_null] = linear_regularization_standardize(X, y, C, L, beta0);

[beta_hat_s, lambda2, ~, SSE, nu, W_s] = do_inversion_standard_aicc (X_s, y_s);

beta_hat = L_generalized_inv*beta_hat_s + (beta_null + beta0);
W = L_generalized_inv*W_s*L_generalized_inv';

mu = Xstar*beta_hat;

[p, n] = size(L);
nu -= n - p; #include unregularized fitted parameters

s2 = SSE / nu;
ss = sum(Xstar.* (W*Xstar')', 2);
sigma = sqrt(s2 * (1 + ss));

endfunction

%{
Example:
m = 100;
n = 6;
p = 5;
b = [100 1 2 3 4 5]';
beta0 = zeros(n, 1);
C = diag(exp((randn(m, 1) + 1) .^ 2)); #highly heteroscedastic noise
L = [zeros(p, n-p) eye(p, p)];
X = randn(m, n);
y = X*b + randn(m, 1) .* sqrt(diag(C));
[X_s y_s L_generalized_inv x_null] = linear_regularization_standardize(X, y, C, L, beta0);
[X_s,y_s,L_p,K,M] = std_form(X,L,y); %function from regtools

[beta_hat_s, lambda2, r, SSE, nu, W_s] = do_inversion_standard_aicc (X_s, y_s);
beta_hat = L_generalized_inv*beta_hat_s + (x_null + beta0);
W = L_generalized_inv*W*L_generalized_inv'; #covariance of original variables
Xstar = randn(1, n);

#these both give the same result
[mu1, sigma1, nu1] = regress_reg_gen_aicc_params (X, y, Xstar, speye(m), speye(n), zeros(n, 1));
[mu2, sigma2, nu2] = regress_reg_aicc_params (X, y, Xstar);
#this should be a bit different (slightly smaller uncertainty because more variability can be explained by regression parameters)
[mu3, sigma3, nu3] = regress_reg_gen_aicc_params (X, y, Xstar, speye(m), L, zeros(n, 1));
#smaller uncertainty because can better account for heteroscedastic noise in training data
[mu4, sigma4, nu4] = regress_reg_gen_aicc_params (X, y, Xstar, C, speye(n), zeros(n, 1));
Xstar*b #"actual" value
%}
