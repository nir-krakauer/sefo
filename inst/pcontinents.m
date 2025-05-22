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

## PCONTINENTS 	Plots continent outlines with longitude from -180-+180 deg.
## See also pcontinents2

## Author: Nir Krakauer <nkrakauer@ccny.cuny.edu>

function [varargout] = pcontinents


  # load percent-land raster (from topo.mat) 
  st = warning ("query", "Octave:data-file-in-path").state;
  warning ("off", "Octave:data-file-in-path")
  load topo
  warning (st, "Octave:data-file-in-path")
  
  
  # plot shoreline
  hold on
  [cs,h] = contour (-179.5:179.5, -89.5:89.5, flipud(topo)-50, [0 0], 'k-');
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
  #xlim([-180 180])
