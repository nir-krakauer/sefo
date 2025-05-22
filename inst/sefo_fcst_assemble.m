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

##read past/current NMME forecast ensemble means, and fill in any missing values 

## Author: Nir Krakauer <nkrakauer@ccny.cuny.edu>



function file = sefo_fcst_assemble (predict_year, predict_month, lag, opts)

  data_dir = getfield (opts, 'data_dir');
  re_download = getfield (opts, 're_download');
  re_process = getfield (opts, 're_process');  
  var_name = getfield (opts, 'var_name');
  fcst_lons = getfield (opts, 'fcst_lons');
  fcst_lats = getfield (opts, 'fcst_lats');
  fcst_start_year = getfield (opts, 'fcst_start_year');
  p_spline = getfield (opts, 'p_spline');
  gcms = getfield (opts, 'gcms');
  
  file = [data_dir 'nmme_' var_name '_' num2str(predict_year) '_' num2str(predict_month, '%2.2i') '_lag_' num2str(lag, '%2.2i') '.mat'];


  if ~exist(file, "file") || re_process

    n_gcms = numel (gcms);
    ny = predict_year - fcst_start_year + 1;
    np = numel(fcst_lats)*numel(fcst_lons);

    
    fcsts_mean = fcsts_sd = nan (np, ny, n_gcms, "single");
    nrs = zeros (ny, n_gcms, "single");
    
    for y = 1:ny
      year = fcst_start_year + y - 1;
      for g = 1:n_gcms
        gcm = gcms{g};
        file_g = sefo_fcst_read (year, predict_month, lag, gcm, opts);
        load(file_g)
        if isempty(vals_mean)
          warning ([file_g ' empty'])
        else
          fcsts_mean(:, y, g) = vals_mean'(:);
          fcsts_sd(:, y, g) = vals_sd'(:);
          nrs(y, g) = nr;
        endif
      endfor
    endfor
    
    #exclude models that have some grid points with no valid values
    missing_models_excluded = any(all(~isfinite(fcsts_mean), 2), 1)(:);
    fcsts_mean = fcsts_mean(:, :, ~missing_models_excluded);
    fcsts_sd = fcsts_sd(:, :, ~missing_models_excluded);
    nrs = nrs(:, ~missing_models_excluded);
    gcms = gcms(~missing_models_excluded);
    n_gcms = numel (gcms);
    
    #fill in missing values
    fcsts_sd_mean = rms (fcsts_sd, 2);
    for y = 1:ny
      missing_models = isnan(squeeze(fcsts_mean(1, y, :)));
      n_missing = sum (missing_models);
      n_nonmissing = sum (~missing_models);
      if any(missing_models)
        if all(missing_models)
          error (['No forecasts found for ' num2str(predict_month, '%2.2i') '/' num2str(predict_year) ' at lag ' num2str(lag, '%2.2i')])
        endif
        warning ([file ', year ' num2str(y) ': filling in missing forecasts for model(s) ' num2str(find(missing_models)(:)')])
        y_avail = all(isfinite(squeeze(fcsts_mean(1, :, :)))');
        if n_nonmissing == 1
          adj_factor = mean(fcsts_mean(:, y_avail, missing_models) - fcsts_mean(:, y_avail, ~missing_models), 2);
          fcsts_mean(:, y, missing_models) = fcsts_mean(:, y, ~missing_models) + adj_factor;
        else
          adj_factor = mean(fcsts_mean(:, y_avail, missing_models) - mean(fcsts_mean(:, y_avail, ~missing_models), 3), 2);
          fcsts_mean(:, y, missing_models) = mean(fcsts_mean(:, y, ~missing_models), 3) + adj_factor;
        endif
        fcsts_sd(:, y, missing_models) = fcsts_sd_mean(:, missing_models);
      endif
    endfor
    
    #add trends
    fcsts_mean_trend = nan (np, ny, n_gcms, "single");
    for m = 1:n_gcms
      fcsts_mean_trend(:, :, m) = csaps((1:ny)', double(fcsts_mean(:, :, m)'), p_spline, (1:ny)');
    endfor
    
    try
      save (file, "fcsts_mean", "fcsts_sd", "fcsts_mean_trend", "nrs", "gcms", "missing_models_excluded")
    catch
      warning(["error saving " file])
      save ("-v6", file, "fcsts_mean", "fcsts_sd", "fcsts_mean_trend", "nrs", "gcms", "missing_models_excluded")
    end_try_catch
  endif
  
endfunction

%{
plot(fcst_start_year:predict_year, squeeze(mean(fcsts_mean, 1))); grid on

[yi, p, sigma2] = csaps_sel((1:ny)', reshape(permute(fcsts_mean, [2 1 3]), ny, np*n_gcms), (1:ny)', [], 'gcv');
#aicc: 3.7377e-09, vm: 3.7377e-09, gcv: 4.6082e-05
plot((1:ny)', center(fcsts_mean_trend(1:1001:end, :), 2)')

z = mean(fcsts_mean(:, end, :) - fcsts_mean_trend(:, end, :), 3);
z = mean(fcsts_mean(:, end, :), 3);
colormap('default'); imagesc(fcst_lons, fcst_lats, reshape(z, numel(fcst_lons), numel(fcst_lats))', [min(z) max(z)]); axis xy; colorbar; pcontinents2
%}

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
%!                 'fcst_start_year', 1982, ...
%!                 'p_spline', 3E-4, ...   
%!                 'gcms', {{'CMC1-CanCM3', 'CMC2-CanCM4', 'COLA-RSMAS-CCSM3', 'COLA-RSMAS-CCSM4', 'GFDL-CM2p1-aer04', 'GFDL-CM2p5-FLOR-A06', 'GFDL-CM2p5-FLOR-B01', 'NASA-GMAO-062012', 'NCEP-CFSv2'}} ...
%!                 );
%! predict_year = 2015; predict_month = 4; lag = 0;
%! file = sefo_fcst_assemble (predict_year, predict_month, lag, opts)
%! load(file, "fcsts_mean"); disp(squeeze(mean(fcsts_mean)))
       
