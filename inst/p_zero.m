## Copyright (C) 2016 Nir Krakauer
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

#probability of mean of (possibly autocorrelated) time series (row vector or matrix of row vectors), under null hypothesis of zero mean
#and the minimum mean difference needed to produce a significance level of alpha
#
#reference: NY Krakauer, MJ Puma, BI Cook (2013), Impacts of soil-aquifer heat and water fluxes on simulated global climate, Hydrology and Earth System Sciences, 17(5): 1963-1974, doi: 10.5194/hess-17-1963-2013

function [p, dt] = p_zero (y, alpha=0.05)
  n = size (y, 2);
  m = mean (y, 2);
  s = std (y, [], 2);
  y = y - m;
  c = sum (y(:, 1:(end-1)) .* y(:, 2:end), 2) / (n-1);
  r = c ./ (s .^ 2);
  n_adj = n * (1 - r) ./ (1 + r);

  if isargout (1)
    t_score = sqrt(n_adj) .* (m ./ s);
    p = 2 * tcdf (-abs(t_score), n_adj-1);
  endif

  if isargout (2)
    dt = (s ./ sqrt(n_adj)) .* tinv ((1 - alpha/2), n_adj-1);
  endif
  
endfunction
