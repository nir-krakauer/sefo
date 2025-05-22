## Copyright (C) 2009-2017 Nir Krakauer
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

%creates a color map with a specified number of colors that goes from blue to red, with gray in the middle
%this is good for dichotomizing the top and bottom halves of the range better
%one input argument: number of desired colors (if none is provided, this is taken to be the number of colors in the current color map, which by default is 64)
%one output: the generated color map

%based on Kai Habel's function jet

## Author: Nir Krakauer <nkrakauer@ccny.cuny.edu>


function map = bluered (number)

  if (nargin == 0)
    number = rows (colormap);
  elseif (nargin == 1)
    if (! isscalar (number))
      error ("bluered: argument must be a scalar");
    endif
  else
    print_usage ();
  endif

  graylevel = 0.9; %set to something between 0 and 1 -- the lower this is, the darker the middle color

if (number > 1)
    x = linspace(0, 1, number)';
    r = g = b = zeros(number, 1);
    inds = x < 0.5 & x > 0.25;
    r(inds) = 2*x(inds) - 0.5;
    inds = x > 0.5;
    r(inds) = 1;
    g(inds) = 2*(1 - x(inds));
    inds = x < 0.5;
    g(inds) = 2*x(inds);
    b(inds) = 1;
    map = [r, g, b];
    inds = (x == 0.5);
    map(inds, :) = graylevel;
else
    map = graylevel * ones(number, 3);  #just gray
endif

endfunction

