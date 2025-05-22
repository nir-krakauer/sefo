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

##  CRPS with a Student t distribtion
##
##  h = t_crps(X,DF);
##
##  Computes the Continuous Ranked Probability Score (CRPS) relative to X
##  of the Student t distribution with DF degrees of freedom.
##
##  For a non-standardized t distribution with mean MU and scale SIGMA,
##  CRPS = SIGMA * t_crps((X-MU)/SIGMA,DF)
##
##  Reference: Alexander Jordan, Closed form expressions for the continuous ranked probability score, https://github.com/FK83/scoringRules/blob/master/crps.pdf

## Author: Nir Krakauer <mail@nirkrakauer.net>

function c = t_crps(x,df)

 c = 2*tpdf(x,df).*(df + x.^2)./(df - 1) + x .* (2*tcdf(x,df)-1) - 2*(sqrt(df)./(df-1)).*exp(betaln(0.5,df-0.5)-2*betaln(0.5,df/2));


%!assert(t_crps([-1.2 0 0.4 1 2 3], 10), [0.7426977 0.2447397 0.3061025 0.6022684 1.4235096 2.3881591], 1E-7)
  #[compare with crps(y = c(-1.2, 0, 0.4, 1, 2, 3), family = "t", location = 0, scale = 1, df = 10) in the R scoringRules package]
