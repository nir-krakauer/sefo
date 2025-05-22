
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

#t-distribution parameters (mu, sigma, nu) for MA(n) fit based on the last n values
#syntax: [mu, sigma, nu] = ma_params (y, n)

## Author: Nir Krakauer <nkrakauer@ccny.cuny.edu>




function [mu, sigma, nu] = ma_params (y, n)

ny = size(y, 1);
y = y(max(1, ny-n+1):end, :);

mu = mean(y);
sigma = std(y);
nu = min(ny, n) - 1;
