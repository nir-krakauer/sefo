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

##Obtain and map temperature anomaly forecast for next month (by default) based on current NMME output.
##Can initialize predict_year, predict_month, lag to get forecasts for other months.

## Author: Nir Krakauer <nkrakauer@ccny.cuny.edu>


if !exist("verbose", "var")
  verbose = true;
endif

if verbose
  page_output_immediately(true);
  page_screen_output(false);
  disp("sefo_example: configuring settings...")
endif

if !exist("generate_more_figures", "var")
  generate_more_figures = false;
endif

[this_year this_month this_day] = datevec (date);

if !exist("predict_year", "var") || !exist("predict_month", "var") || !exist("lag", "var")
  months_from_now = 2;

  if this_day > 8 #current month's batch of NMME forecasts is assumed to be uploaded
    lag = months_from_now;
  else #use previous month's batch
    lag = 1 + months_from_now;
  endif

  predict_month = this_month + months_from_now;
  predict_year = this_year;
  if predict_month > 12
    predict_month -= 12;
    predict_year += 1;
  endif

  #predict_year = 2015; predict_month = 12; lag = 1;
endif

if !exist("data_dir", "var")
  data_dir = [pwd filesep()];
endif

if !exist("plot_dir", "var")
  plot_dir = [pwd filesep()];
endif

opts = struct ('data_dir', data_dir, ...
                're_download', false, ...
                're_process', false, ...
                'var_name', 'tref', ...
                'fcst_lons', (0:359)', ...
                'fcst_lats', (-90:90)', ... 
                'fcst_start_year', 1982, ...
                'obs_start_year', 1957, ...  
                'obs_end_year', min(this_year, predict_year)-1, ...
                'obs_dataset', 'best', ...
                'p_spline', 1E-4, ...
                'dof_spline', 2.5, ...
                'lamda', 10, ...
                'predict_adj', false, ...
                'predict_adj_years', 1994:1998, ...
                'predict_adj_months', 1:12, ...                
                'methods', {{'EW-a', 'MM'}}, ...
                'clim_years', 1981:2010, ...
                'gcms', {{'CMC1-CanCM3', 'CMC2-CanCM4', 'COLA-RSMAS-CCSM4', 'GFDL-CM2p1-aer04', 'GFDL-CM2p5-FLOR-A06', 'GFDL-CM2p5-FLOR-B01', 'NASA-GMAO-062012', 'NCEP-CFSv2'}} ...
                );

if verbose
  disp(["sefo_example: downloading and processing observations and forecasts (this may take some minutes; around 1 GB space will be used in " data_dir ")..."])
  disp("")
endif
                
method = 'EW-a';
file = sefo_predict (predict_year, predict_month, lag, method, opts){1}; load(file);
method_trend = method;
mu_trend = mu;
sigma_trend = sigma;
nu_trend = nu;

method = 'MMT';
file = sefo_predict (predict_year, predict_month, lag, method, opts){1}; load(file);
method_model_trend = method;
mu_model_trend = mu;
sigma_model_trend = sigma;
nu_model_trend = nu;

if verbose
  disp("sefo_example: generating map of calibrated forecast...")
endif

title_string = ['Forecast temperature for ' num2str(predict_month) '/' num2str(predict_year) ' relative to trend, K'];

obs_start_year = getfield (opts, 'obs_start_year');
clim_years = getfield (opts, 'clim_years');
lons = getfield (opts, 'fcst_lons');
lats = getfield (opts, 'fcst_lats');
nlats = numel (lats);
nlons = numel (lons);
obs_dataset = getfield (opts, 'obs_dataset');

file_obs = sefo_obs_assemble (predict_year, predict_month, opts);
load (file_obs, "obs", "obs_inds")


if strcmp (obs_dataset, "bestland") #interpolate to the 1*1 grid
  m = sum (obs_inds);
  n = numel (obs_inds);

  source_file_name = ['Complete_TAVG_EqualArea.nc'];
  source_file = [data_dir source_file_name];
  lon_obs = ncread (source_file, 'longitude');
  lat_obs = ncread (source_file, 'latitude');
  lons_arr = ones(nlats, 1) * lons';
  lats_arr = lats * ones(1, nlons);

  A = sparse (n, m);
  scale = 1;
  max_scale = 8;
  max_dist = 1.5;
  for i = 1:m
    d = sphere_dist (lon_obs(i), lat_obs(i), lons_arr, lats_arr);
    inds = find (d < max_scale);
    w = exp (-d(inds) ./ scale);
    A(inds, i) = double (w);
  endfor
  w_min = exp (-max_dist ./ scale);
  obs_inds = full (max (A, [], 2) >= w_min);
  n = sum (obs_inds);
  A = A(obs_inds, :);
  a = full (sum (A, 2));
  A = spdiags(1 ./ a, 0, n, n) * A;
  #for i = 1:n
  #  A(i, :) ./= a(i); #A ./= sum(A, 2) doesn't work with A sparse
  #endfor
  obs_orig = obs;
  ny = size (obs, 2);
  obs = nan (n, ny);
  for i = 1:ny
    obs(:, i) = A * obs_orig(:, i);
  endfor
  mu_trend = A*mu_trend;
  sigma_trend = A*sigma_trend;
  if !isscalar (nu_trend)
      nu_trend = A*nu_trend;
  endif
  mu_model_trend = A*mu_model_trend;
  sigma_model_trend = A*sigma_model_trend;
  if !isscalar (nu_model_trend)
      nu_model_trend = A*nu_model_trend;
  endif
endif

arr = nan (nlats, nlons);
arr(obs_inds) = mu_model_trend - mu_trend;
ll = max(abs(quantile(arr(obs_inds), [0.05 0.95])));
colormap("default"); imagesc_withnan (lons, lats, arr, [-ll ll]); axis xy; colorbar; pcontinents2
title (title_string)

h = gcf;
FS = findall (h,'-property','FontSize');
set (FS,'FontSize', 16);

plot_name_string = [num2str(predict_year) '_' num2str(predict_month, '%2.2i') '_lag_' num2str(lag, '%i')]; 
plot_file = ['model_trend_difference_' plot_name_string];

if !generate_more_figures
  return
endif

#to generate and save additional maps comparing the trend and model ensemble based forecasts:
print('-dpng', [plot_dir plot_file '.png'])

sd_trend = sigma_trend .* sqrt(nu_trend ./ (nu_trend - 2));
sd_model_trend = sigma_model_trend .* sqrt(nu_model_trend ./ (nu_model_trend - 2));
sd_reduction = 1 - (sd_model_trend ./ sd_trend);

clim_mean = mean (obs(:, clim_years - obs_start_year + 1), 2);
arr(obs_inds) = mu_trend - clim_mean;
ll = max(abs(quantile(arr(obs_inds), [0.05 0.95])));
title_string = ['Trend temperature for ' num2str(predict_month) '/' num2str(predict_year) ' relative to ' num2str(clim_years(1)) '-' num2str(clim_years(end)) ' climatology, K'];
colormap("default"); imagesc_withnan (lons, lats, arr, [-ll ll]); axis xy; colorbar; pcontinents2
title (title_string)
plot_file = ['trend_difference_' plot_name_string];
print('-dpng', [plot_dir plot_file '.png'])

arr(obs_inds) = mu_model_trend - clim_mean;
ll = max(abs(quantile(arr(obs_inds), [0.05 0.95])));
title_string = ['NMME temperature for ' num2str(predict_month) '/' num2str(predict_year) ' relative to ' num2str(clim_years(1)) '-' num2str(clim_years(end)) ' climatology, K'];
colormap("default"); imagesc_withnan (lons, lats, arr, [-ll ll]); axis xy; colorbar; pcontinents2
title (title_string)
plot_file = ['model_difference_' plot_name_string];
print('-dpng', [plot_dir plot_file '.png'])

arr(obs_inds) = 100*sd_reduction;
ll = max(abs(quantile(100*sd_reduction, [0.05 0.95])));
title_string = ['Percent reduction in forecast SD for ' method_model_trend ' versus ' method_trend];
colormap("default"); imagesc_withnan (lons, lats, arr, [-ll ll]); axis xy; colorbar; pcontinents2
title (title_string)
plot_file = ['sd_reduction_' plot_name_string];
print('-dpng', [plot_dir plot_file '.png'])

qs = [1/3 2/3];
nqs = numel (qs);
vals = quantile (obs(:, clim_years - obs_start_year + 1), qs, 2);
cdf_vals = nan (size(vals));
for v = 1:nqs   
  cdf_vals(:, v) = tcdf ((vals(:, v) - mu_trend) ./ sigma_trend, nu_trend);
endfor
arr(obs_inds) = 0;
for d = 6:-1:1
 arr(find(obs_inds)(cdf_vals(:, 1) > (1 - d/10) & (cdf_vals(:, 1) > (1 - cdf_vals(:, 2))))) = d-7;
endfor
for d = 6:-1:1
 arr(find(obs_inds)(cdf_vals(:, 2) < (d/10) & (cdf_vals(:, 1) < (1 - cdf_vals(:, 2))))) = 7-d;
endfor
c = colormap(bluered(13)); c(7, :) = 1; colormap(c); #make white instead of black the middle color
imagesc_withnan(lons, lats, arr, [-6 6]); axis xy; colorbar; pcontinents2
cbh = findobj( gcf(), 'tag', 'colorbar');
set( cbh,                
'linewidth', 2,
'tickdir', 'out',
'ticklength',[0.005,0.005],
'ytick', [-4.5, -2.5, -0.5, 0.5, 2.5, 4.5],
'yticklabel',{'80%','60%','40%','40%','60%', '80%'});
title_string = [num2str(clim_years(1)) '-' num2str(clim_years(end)) ' climatology T tercile probabilities based on trend for ' num2str(predict_month) '/' num2str(predict_year)];
title (title_string)
plot_file = ['trend_terciles_' plot_name_string];
print('-dpng', [plot_dir plot_file '.png'])

for v = 1:nqs   
  cdf_vals(:, v) = tcdf((vals(:, v) - mu_model_trend) ./ sigma_model_trend, nu_model_trend);
endfor
#arr = zeros(nlats, nlons);
arr(obs_inds) = 0;
for d = 6:-1:1
 arr(find(obs_inds)(cdf_vals(:, 1) > (1 - d/10) & (cdf_vals(:, 1) > (1 - cdf_vals(:, 2))))) = d-7;
endfor
for d = 6:-1:1
 arr(find(obs_inds)(cdf_vals(:, 2) < (d/10) & (cdf_vals(:, 1) < (1 - cdf_vals(:, 2))))) = 7-d;
endfor
c = colormap(bluered(13)); c(7, :) = 1; colormap(c); #make white instead of black the middle color
imagesc_withnan(lons, lats, arr, [-6 6]); axis xy; colorbar; pcontinents2
cbh = findobj( gcf(), 'tag', 'colorbar');
set( cbh,                
'linewidth', 2,
'tickdir', 'out',
'ticklength',[0.005,0.005],
'ytick', [-4.5, -2.5, -0.5, 0.5, 2.5, 4.5],
'yticklabel',{'80%','60%','40%','40%','60%', '80%'});
title_string = [num2str(clim_years(1)) '-' num2str(clim_years(end)) ' climatology T tercile probabilities based on NMME for ' num2str(predict_month) '/' num2str(predict_year)];
title (title_string)
plot_file = ['model_terciles_' plot_name_string];
print('-dpng', [plot_dir plot_file '.png'])


for v = 1:nqs   
  vals(:, v) = mu_trend + tinv(qs(v), nu_trend) .* sigma_trend;
  cdf_vals(:, v) = tcdf((vals(:, v) - mu_model_trend) ./ sigma_model_trend, nu_model_trend);
endfor
arr(obs_inds) = 0;
for d = 6:-1:1
 arr(find(obs_inds)(cdf_vals(:, 1) > (1 - d/10) & (cdf_vals(:, 1) > (1 - cdf_vals(:, 2))))) = d-7;
endfor
for d = 6:-1:1
 arr(find(obs_inds)(cdf_vals(:, 2) < (d/10) & (cdf_vals(:, 1) < (1 - cdf_vals(:, 2))))) = 7-d;
endfor
c = colormap(bluered(13)); c(7, :) = 1; colormap(c); #make white instead of black the middle color
imagesc_withnan(lons, lats, arr, [-6 6]); axis xy; colorbar; pcontinents2
cbh = findobj( gcf(), 'tag', 'colorbar');
set( cbh,                
'linewidth', 2,
'tickdir', 'out',
'ticklength',[0.005,0.005],
'ytick', [-4.5, -2.5, -0.5, 0.5, 2.5, 4.5],
'yticklabel',{'80%','60%','40%','40%','60%', '80%'});
title_string = ['Trend T tercile probabilities based on NMME for ' num2str(predict_month) '/' num2str(predict_year)];
title (title_string)
plot_file = ['model_trend_terciles_' plot_name_string];
print('-dpng', [plot_dir plot_file '.png'])



