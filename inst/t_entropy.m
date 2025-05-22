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

% Entropy of Student t distribtion
%
% h = t_entropy(DF);
%
% Computes the entropy of the Student distribution 
%    with DF degrees of freedom. 
%
% see also: TCDF, TINV 

## Author: Nir Krakauer <nkrakauer@ccny.cuny.edu>

function h = t_entropy(n)

h = (n-1)/2 .* (psi((1 + n)/2) - psi(n/2)) + log(n)/2 + betaln(n/2, 1/2);

%!assert(t_entropy(1E3), (log(2*pi)+1)/2, 1E-6) #converges to the standard normal distribution's entropy when DF is large
