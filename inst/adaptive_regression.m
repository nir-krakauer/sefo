
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

#[ystar, sstar, gamma] = adaptive_regression(y, X, R, gamma_in)
#
#Estimate ystar at Xstar = [1 zeros(k-1, 1)]' given data y ~ N(X*beta, sigma^2 * Q), with beta and sigma2 unknown and Q = (1-gamma)*eye + gamma*R, where if gamma is unknown it is estimated by maximum likelihood.
#
#R may be a covariance matrix or a cell array containing its eigendecomposition {P, D}, with P orthonormal and D diagonal and P*D*P' = R
#
#Gamma can be specified as a single value or a range.
#
#References:
#
#Cooley, T. F. & Prescott, E. C. (1973) An adaptive regression model, International Economic Review, 14
#
#Krakauer, N. Y. & Devineni, N. (2015) Up-to-date probabilistic temperature climatologies, Environmental Research Letters, 10:024014

## Author: Nir Krakauer <nkrakauer@ccny.cuny.edu>



function [ystar, sstar, gamma] = adaptive_regression(y, X, R, gamma_in)

if iscell (R)
  P = R{1};
  D = R{2};
else
  [P, D] = eig (R);
endif

if nargin < 4 || isempty(gamma_in)
  gamma_in = [0 1];
endif
  

[T, k] = size (X);
p = k - 1;

#estimate gamma by maximum likelihood
#L = -(T/2)*s^2 - (1/2)*log(det(Q))
#B = inv(X'*Qinv*X)*X*Qinv*y
#r = y - X*B
#s^2 = r'*Qinv*r / T
#now
#Q = (1-gamma)*eye(T) + gamma*R = (1-gamma)*eye(T) + gamma*P*D*P'
#Qinv = P*inv((1-gamma)*eye(T) + gamma*D)*P';
#det(Q) = prod((1-gamma) + gamma*diag(D))
#log(det(Q)) = sum(log((1-gamma) + gamma*diag(D)))

f_cost = @(gamma) cost(gamma, y, X, T, P, D);
if numel(gamma_in) == 2
  gamma = fminbnd (f_cost, gamma_in(1), gamma_in(2));
else
  gamma = gamma_in;
endif
Qinv = P*inv((1-gamma)*eye(T) + gamma*D)*P';
W = pinv(X'*Qinv*X);
B = W*X'*Qinv*y;
ystar = B(1);

#ystar = Xstar*W*X'*Qinv*y;
sigma2 = (1/(T - p)) * y'*(Qinv - Qinv*X*W*X'*Qinv)*y; #or (1/(T - p)) * (y - ystar)'*Qinv*(y - ystar)
#sstar = sqrt(sigma2 * (1 + W(1, 1))); #prediction SE for y*
sstar = sqrt(sigma2 * ((1 - gamma) + W(1, 1))); #prediction SE for y*

endfunction

function J = cost(gamma, y, X, T, P, D)
  Qinv = P*inv((1-gamma)*eye(T) + gamma*D)*P';
  log_det_Q = sum(log((1-gamma) + gamma*diag(D)));
  B = inv(X'*Qinv*X)*X'*Qinv*y;
  r = y - X*B;
  s2T = r'*Qinv*r;
  J = s2T + log_det_Q; #halve to get negative log likelihood
endfunction

