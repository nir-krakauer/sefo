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

##read historical temperatures regridded to NMME

## Author: Nir Krakauer <nkrakauer@ccny.cuny.edu>




function file = sefo_obs_assemble (predict_year, predict_month, opts)

  data_dir = getfield (opts, 'data_dir');
  re_process = getfield (opts, 're_process');  
  var_name = getfield (opts, 'var_name');
  fcst_lons = getfield (opts, 'fcst_lons');
  fcst_lats = getfield (opts, 'fcst_lats');
  obs_start_year = getfield (opts, 'obs_start_year');
  obs_end_year = getfield (opts, 'obs_end_year');  
  obs_dataset = getfield (opts, 'obs_dataset');
  p_spline = getfield (opts, 'p_spline');

  file = [data_dir 'obs_assemble_' var_name '_' obs_dataset '_' num2str(obs_start_year) '_' num2str(predict_year) '_' num2str(predict_month, '%2.2i') '.mat'];
  
  if ~exist(file, "file") || re_process

    obs_end_year = min(predict_year-1, obs_end_year);
    
    ny = obs_end_year - obs_start_year + 1;
    ny_predict = predict_year - obs_start_year + 1;

    np = numel(fcst_lats)*numel(fcst_lons);  
    obs = nan (np, ny); 
    
    for y = 1:ny
      year = obs_start_year + y - 1;
      file_year = sefo_obs_read (year, predict_month, opts);
      load (file_year)
      obs(:, y) = arr_regrid(:);
    endfor
    
    #keep only the grid points with observations for all years
    obs_inds = all(isfinite(obs), 2);
    
    obs = obs(obs_inds, :); #17343 points for CAMS
    
    obs_trend = csaps ((1:ny)', obs', p_spline, [(1:ny) ny_predict]');
    
    save (file, "obs", "obs_trend", "obs_inds")
  
  endif
  

  
endfunction
       
#colormap('default'); imagesc_withnan(fcst_lons, fcst_lats, reshape(mean(obs, 2), [numel(fcst_lats), numel(fcst_lons)])); axis xy; colorbar
#plot(obs_start_year + (0:(ny-1)), mean(obs))
#[yi, p, sigma2] = csaps_sel((1:ny)', obs', (1:ny)', [], 'aicc');
#plot((1:ny)', center(yi(1:154:end, :), 2)')
#gcv: 7.7284e-04; aicc: 3.9722e-04; vm: 1.2215e-05
#[yi, p, sigma2] = csaps_sel((1:ny)', obs', (1:ny)', [], 'gcv');
#plot((1:ny)', center(obs_trend(1:154:end, :), 2)')


%!demo
%! if !exist("data_dir", "var")
%!   data_dir = [pwd '/'];
%! endif
%! opts = struct ('data_dir', data_dir, ...
%!                 're_download', false, ...
%!                 're_process', false, ...
%!                 'var_name', 'tref', ...
%!                 'fcst_lons', (0:359)', ...
%!                 'fcst_lats', (-90:90)', ... 
%!                 'obs_dataset', 'best', ...
%!                 'obs_start_year', 1957, ...
%!                 'obs_end_year', 2015, ...
%!                 'p_spline', 3E-4 ...   
%!                 );
%! predict_year = 2016; predict_month = 1; lag = 1;
%! file = sefo_obs_assemble (predict_year, predict_month, opts)
