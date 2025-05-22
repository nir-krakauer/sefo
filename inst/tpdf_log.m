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

% TPDF_LOG: logarithm of the student probability density 
%
% pdf = tpdf_log(x,DF);
%
% Computes the logarithm of the PDF of a student t distribution
%    with DF degrees of freedom
% x,DF must be matrices of same size, or any one can be a scalar. 
%
% see also: TPDF, TINV, TCDF, NORMPDF, NORMCDF, NORMINV 

%based on tpdf in the NaN toolbox

## Author: Nir Krakauer <nkrakauer@ccny.cuny.edu>

function p = tpdf_log(x,n)

% allocate memory and check size of arguments
p = x+n;	  % if this line causes an error, size of input arguments do not fit.

% make size of x and n equal
n = x+n-x;
x = x+n-n;

% workaround for invalid input values
ix = (n>0) & (n~=inf) & ~isnan(x);
if any(ix(:))
  p(ix) = -(n(ix)+1).*log(1+x(ix).^2./n(ix))/2 - log(n(ix))/2 - betaln(n(ix)/2, 1/2);
end; 
p(~ix)= NaN;

% shape output
p = reshape (p,size(x));

%!assert(tpdf_log(NaN,4),NaN)
