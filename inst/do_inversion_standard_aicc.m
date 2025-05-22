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

%solves Ax \approx b 
%assuming that each element of x has an independent prior distribution that is normal with zero mean and unknown standard deviation s/lambda, and that the error in each element of b is iid standard normal with unknown standard deviation s
%use the Akaike Information Criterion with second-order correction to select lambda (unless its square is given)
%Syntax: [x, lambda2, r, rss, T, W] = do_inversion_standard_aicc(A, b[, lambda2])

## Author: Nir Krakauer <nkrakauer@ccny.cuny.edu>


function [x, lambda2, r, rss, T, W] = do_inversion_standard_aicc(A, b, lambda2)



    [n, p] = size(A);

    [U, S, V] = svd(A, 0);
    s = diag(S); %singular values
    s2 = s.^2;

    if nargin < 3 || isempty(lambda2)
      aiccf = @(log_lambda2) aicc_function(log_lambda2, A, b, U, s, V, s2, n);
      [log_lambda2 aicc_min] = fminbnd(aiccf, log(s(end)/1E3 + eps), log(s(1)*1E3 + eps));
      lambda2 = exp(log_lambda2);
    endif
    
    x = V*diag(s ./ (s2 + lambda2))*U'*b;
    r = A*x - b; %residuals
    rss = sumsq(r); %residual sum of squares
    T = n - sum(s2 ./ (s2 + lambda2)); %retained degrees of freedom
    
    W = V*diag((s ./ (s2 + lambda2)) .^ 2)*V'; %proportional to covariance matrix of x for given lambda2 [for inverse, reciprocate the singular values]

endfunction
    
function aicc = aicc_function(log_lambda2, A, b, U, s, V, s2, n)

lambda2 = exp(log_lambda2);

x = V*diag(s ./ (s2 + lambda2))*U'*b;
rss = sumsq(A*x - b);

#rss = sumsq((A*V*diag(s ./ (s.^2 + lambda2))*U' - speye(n))*b); #avoid, since A*pinv(A) can be too large

k = sum(s2 ./ (s2 + lambda2)); %equivalent degrees of freedom of the fit

aicc = sum(n*(log(rss/n) + 2*k/(n-k-1)));

endfunction

%{
Example:
n = 2000;
k = 20;
lambda_in = 10;
x_in = randn(k, 1) / lambda_in;
s_in = 0.05;
A = randn(n, k);
b = A * x_in + s_in*randn(n, 1);
[x, lambda2, r, rss, T, W] = do_inversion_standard_aicc(A, b);
plot(x_in, x - x_in, '+')
%}
