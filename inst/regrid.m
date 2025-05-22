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

#regridding in d dimensions
#
#syntax: 
#v = regrid(coords, u, new_coords)
#
#coords: cell array containing cell limit coordinate vectors along each of d dimensions (in order)
#u: array of values at the cells to regrid
#new_coords: cell array containing new cell limit coordinate vectors along each of d dimensions (in order)

## Author: Nir Krakauer <nkrakauer@ccny.cuny.edu>


function v = regrid (coords, u, new_coords)

nd = numel (coords);

v = u;

 warn_state = warning ("query", "Octave:broadcast").state;
 warning ("off", "Octave:broadcast"); #turn off warning message for automatic broadcasting
 unwind_protect

for d = 1:nd
    c = coords{d}(:);
    nc = new_coords{d}(:);

    da = max(min(nc(2:end)', c(2:end)) - max(nc(1:(end-1))', c(1:(end-1))), 0);

    da = da ./ sum(da);

    sv = size(v);
    sn = numel(nc) - 1;

    v = ndmult (da, v, [1 d]);

    v = reshape (v, [sn sv(1:d-1) sv(d+1:end)]); #this makes sure to keep leading singleton dimensions in output, which disappear when ndmult is applied
endfor

 unwind_protect_cleanup
 warning (warn_state, "Octave:broadcast");
 end_unwind_protect

if nd > 1
  v = permute(v, [nd:-1:1]);
endif

%{


sv = nan(nd, 1); #size of the output v
for d = 1:nd
  sv(d) = numel(new_coords{d}(:)) - 1;
endfor

u, w: m*n (for example -- nd = 2)
coords: (m+1), (n+1)
new_coords: (k+1), (l+1)
v: k*l

#for regridding with weights:
v = regrid(coords, u .* w, new_coords);
vw = regrid(coords, w, new_coords);
vwr = v ./ vw;


%}



%!demo
%! lons = 0:30:360;
%! lats = -90:18:90;
%! coords = {lons, lats'};
%! new_coords = {0:40:360, -90:30:90};
%! u = (cosd((15:30:345)(:)) * sind((-81:18:81)(:))');
%! v2 = (cosd((20:40:340)(:)) * sind((-75:30:75)(:))');
%! vr = regrid(coords, u, new_coords);

%!demo
%! lons = 0:30:360;
%! lats = -90:18:90;
%! z = 0:50:200;
%! coords = {lons, lats, z};
%! new_coords = {0:40:360, -90:30:90, 0:100:200};
%! u = repmat(cosd((15:30:345)(:)) * sind((-81:18:81)(:))', [1 1 4]); 
%! v2 = repmat(cosd((20:40:340)(:)) * sind((-75:30:75)(:))', [1 1 2]);
%! vr = regrid(coords, u, new_coords);

%!demo
%! lons = 0:30:360;
%! lats = -90:18:90;
%! z = 0:50:200;
%! coords = {lons, lats, z};
%! new_coords = {[0 360], -90:30:90, 0:100:200};
%! u = repmat(cosd((15:30:345)(:)) * sind((-81:18:81)(:))', [1 1 4]); 
%! v2 = repmat(cosd((180)(:)) * sind((-75:30:75)(:))', [1 1 2]);
%! vr = regrid(coords, u, new_coords);