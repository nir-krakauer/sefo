## Copyright (C) 2015-2016 Nir Krakauer
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

## PCONTINENTS2 	Plots continent outlines with longitude from 0-360 deg.
## See also pcontinents

## Author: Nir Krakauer <nkrakauer@ccny.cuny.edu>

function [varargout] = pcontinents2

  # load percent-land raster (from topo.mat)
  st = warning ("query", "Octave:data-file-in-path").state;
  warning ("off", "Octave:data-file-in-path")
  load topo
  warning (st, "Octave:data-file-in-path")
  
  # plot shoreline
  hold on
  [cs,h] = contour (0.5:359.5, -89.5:89.5, topo-50, [0 0], 'k-');
  set (h, 'LineWidth', 0.8);
  grid off
  hold off
  
  if nargout == 1
    varargout{1} = h;
  elseif nargout == 2
    varargout{1} = cs;
    varargout{2} = h;
  endif
  
  #ylim([-90 90])
  #xlim([0 360])
  
  
%{
generating topo.mat:
download ftp://daac.ornl.gov/data/islscp_ii/ancillary/combined_ancillary_xdeg/data/land_water_masks_xdeg.zip
extract land_percent2_1d.asc
topo = dlmread ("land_percent2_1d.asc", " ", 6, 0);
save topo.mat topo
%}
