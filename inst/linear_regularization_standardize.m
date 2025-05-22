
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

# Transform a linear inverse problem with a linear regularization term to standard form (Hansen 1998, Krakauer et al. 2004), i.e. go from
# min_x (Ax - b)'*C^{-1}*(Ax - b) + (L*(x - x0))'*(L*(x - x0))
# to
# min_{x_s} (A_s x_s - b_s)'*(A_s x_s - b_s) + (x_s)'*(x_s)
#
#TODO: can it be generalized to when there are multiple regularization parameters that apply to different subsets of x?

## Author: Nir Krakauer <nkrakauer@ccny.cuny.edu>





function [A_s b_s L_generalized_inv x_null] = linear_regularization_standardize(A, b, C, L, x0)

b -= A*x0; %convert to the equivalent problem with x0 = 0 

n = size(A, 2);
N = pinv (A*(speye(n) - L \ L), (norm(A) + norm(L))/1E10);
L_generalized_inv = (speye(n) - N * A) / L; %A-weighted generalized inverse of L, Hansen Eq. 2.32
x_null = N * b; %Component of x in the null space of L, Hansen Eq. 2.33
A_s = A*L_generalized_inv; %Hansen Eq. 2.35
b_s = b - A*x_null;
b_se = chol(C); %C^{1/2}
A_s = b_se \ A_s; %Hansen Eq. 5.2
b_s = b_se \ b_s;

%Note: x = L_generalized_inv*x_s + x_null + x0; %Cf. Hansen Eq. 2.36

endfunction

%{
Example:
m = 100;
n = 6;
x = [-100 0.01 0.3 0.2 0.1 -0.1]';
A = zeros(m, n);
A(:, 1) = 1;
A(:, 2) = (1:m)';
A(:, 3:6) = 200 + randn(m, 4);
b = A*x + randn(m, 1)/10;
L = [zeros(4, 2) eye(4)];
C = speye(m);
beta0 = zeros(n, 1);
[A_s b_s L_generalized_inv x_null] = linear_regularization_standardize (A, b, C, L, beta0);
[beta_hat_s, ~, ~, SSE, nu, W_s] = do_inversion_standard_aicc (A_s,b_s);
beta_hat = L_generalized_inv*beta_hat_s + (x_null + beta0);
regression_eff = nsc(A*x, A*beta_hat) #0.998 or so

%}
