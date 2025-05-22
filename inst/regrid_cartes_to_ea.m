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

##regrid from a Cartesian to an equal-area latitude-longitude grid

## Author: Nir Krakauer <nkrakauer@ccny.cuny.edu>






function arr_out = regrid_cartes_to_ea(arr, lons_cartes, lats_cartes_lims, lats_ea_lims, lons_ea, ea_lat_inds)
  #first change latitudes
  [nlat_in nlon_in] = size(arr);
  nlat_out = numel(lats_ea_lims) - 1;
  arr_new = nan(nlat_out, nlon_in);
  for j = 1:nlon_in
    arr_new(:, j) = regrid({lats_cartes_lims}, arr(:, j), {lats_ea_lims});
  endfor
  
  #now regrid over longitudes at each new latitude
  arr_out = nan(size(lons_ea));
  lons_cartes = lons_cartes(:)';
  tmp = [lons_cartes-360 lons_cartes lons_cartes+360]; 
  lons_cartes_lims = (tmp(1:(end-1))+tmp(2:end))/2;
  for i = 1:nlat_out
    tmp = lons_ea(ea_lat_inds == i)(:)';
    tmp = [tmp(end)-360 tmp tmp(1)+360];
    lons_ea_lims = (tmp(1:(end-1))+tmp(2:end))/2;
    tmp = [arr_new(i, :) arr_new(i, :) arr_new(i, :)];
    tmp = tmp(2:(end-1))';
    arr_out(ea_lat_inds == i) = regrid({lons_cartes_lims}, tmp, {lons_ea_lims});
  endfor

endfunction
