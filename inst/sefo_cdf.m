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

##returns predicted probability of not exceeding specified values or climatology quantiles

## Author: Nir Krakauer <nkrakauer@ccny.cuny.edu>



function file = sefo_cdf (predict_year, predict_month, lag, method, vals='tercile_map', convert_quantile=false, opts)

  data_dir = getfield (opts, 'data_dir');
  re_download = getfield (opts, 're_download');
  re_process = getfield (opts, 're_process');  
  var_name = getfield (opts, 'var_name');
  fcst_lons = getfield (opts, 'fcst_lons');
  fcst_lats = getfield (opts, 'fcst_lats');
  fcst_start_year = getfield (opts, 'fcst_start_year');
  obs_start_year = getfield (opts, 'obs_start_year');
  obs_dataset = getfield (opts, 'obs_dataset');
  clim_years = getfield (opts, 'clim_years');
  predict_adj = getfield (opts, 'predict_adj');
  predict_adj_years = getfield (opts, 'predict_adj_years');
  predict_adj_months = getfield (opts, 'predict_adj_months');  

  if predict_adj
    adj_tag = '_adj';
  else
    adj_tag = '';
  endif
  
  file = [data_dir 'cdf_' var_name '_' method '_' obs_dataset '_' num2str(predict_year) '_' num2str(predict_month, '%2.2i') '_lag_' num2str(lag, '%2.2i') adj_tag '.mat'];

  if ~exist(file, "file") || re_process

    file_obs = sefo_obs_assemble (predict_year, predict_month, opts);
    load (file_obs)
    file_predict = sefo_predict (predict_year, predict_month, lag, method, opts){1};
    load (file_predict, "mu", "sigma", "nu")
    
    #areas represented by each grid cell
    nlats = numel(fcst_lats);
    nlons = numel(fcst_lons);
    fcst_lat_lims = nan(nlats+1, 1);
    fcst_lat_lims(2:(end-1)) = (fcst_lats(1:(end-1)) + fcst_lats(2:end))/2;
    fcst_lat_lims(1) = -90;
    fcst_lat_lims(end) = 90;
    fcst_lat_lims = sind (fcst_lat_lims);
    area_wgts = diff(fcst_lat_lims) * ones(1, nlons);
    
    #only consider grid cells with available historical observations
    if exist("obs_inds", "var") && ~isempty(obs_inds)
      area_wgts = area_wgts(obs_inds);
    else
      obs_inds = [];
    endif
        
    area_wgts = area_wgts(:);
    area_wgts /= sum(area_wgts);    
    np = numel (mu);
    
    if strcmpi (vals, 'tercile_map') #map the tercile probabilities (red shades for higher-than-equal-chances of being in the top tercile of climatology, blue for higher-than-equal-chances of being in the bottom tercile)

      vals = quantile (obs(:, clim_years - obs_start_year + 1), [1/3 2/3], 2);
      for v = 1:2   
        cdf_vals(:, v) = tcdf((vals(:, v) - mu) ./ sigma, nu);
      endfor
      arr = zeros(nlats, nlons);
      for d = 6:-1:1
       arr(cdf_vals(:, 1) > (1 - d/10) & (cdf_vals(:, 1) > (1 - cdf_vals(:, 2)))) = d-7;
      endfor
      for d = 6:-1:1
       arr(cdf_vals(:, 2) < (d/10) & (cdf_vals(:, 1) < (1 - cdf_vals(:, 2)))) = 7-d;
      endfor
      c = colormap(bluered(13)); c(7, :) = 1; colormap(c); #make white instead of black the middle color
      imagesc_withnan(fcst_lons, fcst_lats, arr, [-6 6]); axis xy; colorbar; pcontinents2
      cbh = findobj( gcf(), 'tag', 'colorbar');
      set( cbh,                
      'linewidth', 2,
      'tickdir', 'out',
      'ticklength',[0.005,0.005],
      'ytick', [-4.5, -2.5, -0.5, 0.5, 2.5, 4.5],
      'yticklabel',{'80%','60%','40%','40%','60%', '80%'});
      plot_file = file(1:(end-3));
      print('-dpng', [plot_file 'png'])
     
    else
    
      nv = size(vals, 2);

      if convert_quantile #convert from climatology quantile to actual values
          qs = vals(:);
          vals = quantile (obs(:, clim_years - obs_start_year + 1), qs, 2);
      endif

      if predict_adj
        adj_file = sefo_adj (predict_adj_years, predict_adj_months, lag, method, opts);
        load(adj_file)
        factors_adj = exp(facs_optim);
        sigma = factors_adj(1) * sigma;
        nu = factors_adj(2) * nu;
      endif 

      cdf_vals = nan (np, nv);
      
      for v = 1:nv   
        cdf_vals(:, v) = tcdf((vals(:, v) - mu) ./ sigma, nu);
      endfor
    
      arr = [];
    
    endif
    
    save (file, "cdf_vals", "arr")
    
  endif
  
  
  
endfunction

%{
fcst_lons = getfield (opts, 'fcst_lons');
fcst_lats = getfield (opts, 'fcst_lats');
nlats = numel(fcst_lats); nlons = numel(fcst_lons);
v = 3; colormap('default'); imagesc(fcst_lons, fcst_lats, 1 - reshape(cdf_vals(:, v), nlats, nlons), [0.25 0.75]); axis xy; colorbar; pcontinents2
mean(cdf_vals)
#cf. with http://www.ncdc.noaa.gov/sotc/service/global/map-percentile-mntp/201502.gif
#and http://www.cpc.ncep.noaa.gov/products/NMME/archive/2015010800/current/tmp2m_Lead1.html
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
%!                 'obs_start_year', 1957, ...  
%!                 'obs_end_year', 2015, ...                  
%!                 'obs_dataset', 'best', ...              
%!                 'p_spline', 3E-4, ...
%!                 'dof_spline', 2.5, ...
%!                'lamda', 10, ...
%!                 'predict_adj', true, ...
%!                 'predict_adj_years', 1994:1998, ...
%!                 'predict_adj_months', 1:12, ...                
%!                 'methods', {{'model_linear_reg', 'multimodel_linear_reg'}}, ...
%!                 'clim_years', 1981:2010, ...   
%!                 'gcms', {{'CMC1-CanCM3', 'CMC2-CanCM4', 'GFDL-CM2p1-aer04', 'COLA-RSMAS-CCSM3', 'NASA-GMAO-062012'}} ...
%!                 );
%! predict_year = 2016; predict_month = 2; lag = 1; method = 'rawmodel'; vals = [0 1/3 1/2 2/3 1]; convert_quantile = true;
%! file = sefo_cdf (predict_year, predict_month, lag, method, 'tercile_map', convert_quantile, opts)       
