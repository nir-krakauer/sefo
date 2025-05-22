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

##regrid from an equal-area to a Cartesian latitude-longitude grid

## Author: Nir Krakauer <nkrakauer@ccny.cuny.edu>



function arr_out = regrid_ea_to_cartes (arr, lats_ea_lims, lons_ea, ea_lat_inds, lons_cartes, lats_cartes_lims)

  nlat_in = numel(lats_ea_lims) - 1;
  nlon_out = numel(lons_cartes);
  nlat_out = numel(lats_cartes_lims) - 1;

  #first regrid over longitudes at each latitude
  arr_new = nan(nlat_in, nlon_out);
  lons_cartes = lons_cartes(:)';
  tmp = [lons_cartes(end)-360 lons_cartes lons_cartes(1)+360]; 
  lons_cartes_lims = (tmp(1:(end-1))+tmp(2:end))/2;
  for i = 1:nlat_in
    tmp = lons_ea(ea_lat_inds == i)(:)';
    tmp = [tmp-360 tmp tmp+360];
    lons_ea_lims = (tmp(1:(end-1))+tmp(2:end))/2;
    tmp = arr(ea_lat_inds == i)(:)';
    tmp = [tmp tmp tmp];
    tmp = tmp(2:(end-1))';
    arr_new(i, :) = regrid({lons_ea_lims}, tmp, {lons_cartes_lims});
  endfor
  
  arr_out = nan(nlat_out, nlon_out);
  for j = 1:nlon_out
    arr_out(:, j) = regrid({lats_ea_lims}, arr_new(:, j), {lats_cartes_lims});
  endfor

endfunction
