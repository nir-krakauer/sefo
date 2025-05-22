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

#t-distribution parameters (mu, sigma, nu) for EW(gamma) fit at xpred
#(with optional polynomial time trend or other covariates, offset to all be 0 at the prediction point)
#gamma may be a vector, in which case vectors are returned for mu and sigma
#syntax: [mu, sigma, nu] = ew_params (x, y, xpred, gamma, d = 0, R = [])

## Author: Nir Krakauer <nkrakauer@ccny.cuny.edu>



function [mu, sigma, nu] = ew_params (x, y, xpred, gamma, d = 0, R = [])

warning ("off", "Octave:broadcast", "local");

n = numel(y);

if isscalar(d)
  X = ones(n, d+1);
  for degree = 1:d
    X(:, degree+1) = (xpred - x(:)) .^ degree;
  endfor
else
  X = d;
  d = size(X, 2);
endif

if isempty(R)
  v = xpred - x(:);
  R = min(v, v');
  [P, D] = eig (R);
endif

mu = sigma = nan (size(gamma));
for i = 1:numel(gamma)
  [ystar, sstar] = adaptive_regression(y, X, R, gamma(i));
  mu(i) = ystar;
  sigma(i) = sstar;
endfor

nu = n - 1 - d;
