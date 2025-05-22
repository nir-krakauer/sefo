## Copyright (C) 2017 Nir Krakauer
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

##Quantify likelihood of setting a new heat record each year, and compare skill of 
##trend vs. NMME -based probabilistic forecasts of this

## Author: Nir Krakauer <nkrakauer@ccny.cuny.edu>

if 0 #generate and score forecasts for 2012-2016

  pkg load netcdf nan

  if !exist("verbose", "var")
    verbose = true;
  endif

  if verbose
    page_output_immediately(true);
    page_screen_output(false);
  endif

  nlons = 360;
  fcst_lons = (0:359)';
  fcst_lats = (-90:90)';
  nlats = numel (fcst_lats);
  fcst_lat_lims = nan(nlats+1, 1);
  fcst_lat_lims(2:(end-1)) = (fcst_lats(1:(end-1)) + fcst_lats(2:end))/2;
  fcst_lat_lims(1) = -90;
  fcst_lat_lims(end) = 90;
  fcst_lat_lims = sind (fcst_lat_lims);
  area_wgts = diff(fcst_lat_lims) * ones(1, nlons);
  area_wgts /= sum(area_wgts(:));

  predict_years = 2012:2016;
  ny_predict = numel (predict_years);
  nm_predict = 12 * ny_predict;


  obs_file_name = "Land_and_Ocean_EqualArea.nc";
  data_url = ['http://berkeleyearth.lbl.gov/auto/Global/Gridded/' obs_file_name];
  urlwrite(data_url, [data_dir obs_file_name]);
  lon = ncread ([data_dir obs_file_name], "longitude");
  lat = ncread ([data_dir obs_file_name], "latitude");
  land = ncread ([data_dir obs_file_name], "land_mask");
  T = ncread ([data_dir obs_file_name], "temperature"); #only 7% missing values
  start_year = 1850;
  nm = size (T, 2);
  np = size (T, 1);

  #find closest NMME grid point to each equal-area BEST point
  nmme_inds = zeros (np, 1, "single");
  X = ones(numel(fcst_lats), 1) * fcst_lons';
  Y = fcst_lats * ones(1, numel(fcst_lons));
  X = X(:); Y = Y(:);
  for p = 1:np
    d = sphere_dist (lon(p), lat(p), X, Y);
    [~, b] = min(d);
    nmme_inds(p) = b;
  endfor
  land_inds = (land >= 0.5);
  na_inds = (lat >= 30) & (lat <= 50) & (lon >= -130) & (lon <= -60) & (land >= 0.5);


  predict_year = predict_years(1);
  predict_month = 0;
  for m = 1:nm_predict

    lag = 0;
    predict_month++;
    if predict_month > 12
      predict_month -= 12;
      predict_year++;
    endif
    m_index = 12*(predict_year - start_year) + predict_month;

    opts = struct ('data_dir', data_dir, ...
                  're_download', false, ...
                  're_process', false, ...
                  'var_name', 'tref', ...
                  'fcst_lons', fcst_lons, ...
                  'fcst_lats', fcst_lats, ...
                  'fcst_start_year', 1982, ...
                  'obs_start_year', 1957, ...  
                  'obs_end_year', max(predict_years), ...
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

    method = 'EW-a';
    file = sefo_predict (predict_year, predict_month, lag, method, opts){1}; load(file);
    method_trend = method;

    #switch to the BEST equal-area grid to simplify subsetting and averaging
    mu = mu(nmme_inds); sigma = sigma(nmme_inds);
    if !isscalar(nu)
      nu = nu(nmme_inds);
    endif

    mu_trend = mu;
    sigma_trend = sigma;
    nu_trend = nu;

    m_prev = (m_index-12):(-12):1;
    record_level = max (T(:, m_prev), [], 2);

    is_record = T(:, m_index) > record_level;

    notrend_prob_record = 1 / (numel(m_prev) + 1);
    trend_prob_record = tcdf ((mu_trend - record_level) ./ sigma_trend, nu_trend);


    nll_notrend = -(log(notrend_prob_record) * is_record + log(1 - notrend_prob_record) * !is_record);
    nll_trend = -(log(trend_prob_record) .* is_record + log(1 - trend_prob_record) .* !is_record);
    nll_mean_trend(m, 1) = mean(nll_trend(land_inds));
    nll_mean_notrend(m, 1) = mean(nll_notrend(land_inds));
    f_notrend(m) = notrend_prob_record; #proportion of record-breaking area expected if climate is stationary
    f_actual(m, 1) = mean(is_record(land_inds)); #actual proportion of record-breaking area
    nll_mean_trend(m, 2) = mean(nll_trend(na_inds));
    nll_mean_notrend(m, 2) = mean(nll_notrend(na_inds));
    f_actual(m, 2) = mean(is_record(na_inds));

    for lag = 0:6
      method = 'MMT';
      file = sefo_predict (predict_year, predict_month, lag, method, opts){1}; load(file);
      method_model_trend = method;
      mu = mu(nmme_inds); sigma = sigma(nmme_inds);
      if !isscalar(nu)
        nu = nu(nmme_inds);
      endif
      mu_model_trend = mu;
      sigma_model_trend = sigma;
      nu_model_trend = nu;



      model_prob_record = tcdf ((mu_model_trend - record_level) ./ sigma_model_trend, nu_model_trend);
      modelp_prob_record = min(max(model_prob_record, min(trend_prob_record)), max(trend_prob_record)); #post-process (truncate) the forecast extremes -- improves performance slightly for a test month, but not for 60 months [at lag 0], so don't use
      nll_model = -(log(model_prob_record) .* is_record + log(1 - model_prob_record) .* !is_record);
      nll_modelp = -(log(modelp_prob_record) .* is_record + log(1 - modelp_prob_record) .* !is_record);
      #note: NLL for a perfect forecast would be zero, so could express IG as relative to a perfect forecast

      nll_mean_model(m, lag+1, 1) = mean(nll_model(land_inds));
      nll_mean_modelp(m, lag+1, 1) = mean(nll_modelp(land_inds));
      nll_mean_model(m, lag+1, 2) = mean(nll_model(na_inds));
      nll_mean_modelp(m, lag+1, 2) = mean(nll_modelp(na_inds));
    endfor
  endfor

  mean(f_notrend)
  mean(f_actual)
  mean(nll_mean_notrend)
  mean(nll_mean_trend)
  mean(nll_mean_model, 1)
  mean(nll_mean_modelp, 1)

endif

%{
#land/NA
0.0060611
0.050294   0.034367
0.26274   0.18136
0.18722   0.15764
0.15890   0.18053   0.18078   0.18216   0.18248   0.18308   0.18512
0.12458   0.14461   0.14645   0.14933   0.14844   0.15356   0.15117
0.15556   0.17546   0.17648   0.17774   0.17841   0.17990   0.18200
0.12448   0.14457   0.14646   0.14929   0.14845   0.15358   0.15118

#global
0.017555
0.061812
mean NLL:
mean(nll_mean_notrend)
0.26714
mean(nll_mean_trend)
0.21408
mean(nll_mean_model, 1) #lag 0:6 -- stays better than trend, presumably because of ability to catch dynamical causes of warm episodes, such as El Ninno
0.15651   0.17899   0.18712   0.19055   0.19308   0.19553   0.19793
mean(nll_mean_modelp, 1) #lag 0:6
0.15809   0.17921   0.18688   0.19007   0.19251   0.19514   0.19743


Plots:
BEST record area fraction

Global normalized skill (IG) for records: trend (horiz line), NMME as function of lag
Same, but for NAm land

%}

if 1 #calculate record area fractions

  #BEST record area fraction
  obs_file_name = "Land_and_Ocean_EqualArea.nc";
  data_dir = "/Users/nirkrakauer/Documents/stuff/data/best/";
  if 0 #to re-download
    data_url = ['http://berkeleyearth.lbl.gov/auto/Global/Gridded/' obs_file_name];
    urlwrite(data_url, [data_dir obs_file_name]);
  endif
  pkg load netcdf
  lon = ncread ([data_dir obs_file_name], "longitude");
  lat = ncread ([data_dir obs_file_name], "latitude");
  land = ncread ([data_dir obs_file_name], "land_mask");
  T = ncread ([data_dir obs_file_name], "temperature"); #only 7% missing values
  start_year = 1850;
  nm = size (T, 2);
  inds = (land >= 0.5);
  frac_record = frac_record_stat = ones (nm, 1);
  for m = 13:nm
    m_prev = (m-12):(-12):1;
    nm_prev = numel (m_prev);
    is_record = (T(inds, m) > max(T(inds, m_prev), [], 2));
    frac_record(m) = sum(is_record) ./ sum(isfinite(T(inds, m)));
    frac_record_stat(m) = 1 / (1 + nm_prev .* mean(isfinite(T(inds, m_prev)(:))));
  endfor
  #for N America
  inds = (lat >= 30) & (lat <= 50) & (lon >= -130) & (lon <= -60) & (land >= 0.5);
  frac_record_na = frac_record_stat_na = ones (nm, 1);
  for m = 13:nm
    m_prev = (m-12):(-12):1;
    nm_prev = numel (m_prev);
    is_record = (T(inds, m) > max(T(inds, m_prev), [], 2));
    frac_record_na(m) = sum(is_record) ./ sum(isfinite(T(inds, m)));
    frac_record_stat_na(m) = 1 / (1 + nm_prev .* mean(isfinite(T(inds, m_prev)(:))));
  endfor

endif

plot_dir = '/Users/nirkrakauer/Documents/stuff/plots/nmme/';
plot_format = 'epsc';
if strcmp (plot_format, 'epsc')
  plot_suffix = 'eps';
else
  plot_suffix = plot_format;
endif
  
lw = 2; %line width
fs = 17; %font size
set (0,'defaulttextfontsize',fs);
set (0,'defaultaxesfontsize',fs)
texx = -0.13; texy = 1.03;

clf('reset');
letter = 'a';
p = semilogy(start_year + ((1:nm)'-0.5)/12, [mavg(frac_record_stat, 60) mavg(frac_record, 60)]);
axis([1930 2018 3E-3 1E-1])
title("Global land")
xlabel('Year')
ylabel('Fraction of records')
set (gca, 'fontsize', fs);
set (p, 'linewidth', lw);
set (gca,'position',[0.15 0.2 0.8 0.7])
l = legend ('Expected', 'Actual', "location", "southwest");
set (l, 'fontsize', fs);
text (texx, texy, ["(" letter ")"], "fontsize", fs, "units", "normalized");
print(['-d' plot_format], [plot_dir 'gl_rec' '.' plot_suffix])

clf('reset');
letter = 'b';
p = semilogy(start_year + ((1:nm)'-0.5)/12, [mavg(frac_record_stat_na, 60) mavg(frac_record_na, 60)]);
axis([1930 2018 3E-3 1E-1])
title("North America")
xlabel('Year')
#ylabel('Fraction of records')
set (gca, 'fontsize', fs);
set (gca, 'linewidth', lw);
text (texx, texy, ["(" letter ")"], "fontsize", fs, "units", "normalized")
set (gca,'position',[0.15 0.2 0.8 0.7]);
print(['-d' plot_format], [plot_dir 'na_rec' '.' plot_suffix])

#land
clf('reset');
letter = 'a';
skill_trend = 1 - 0.18722/0.26274;
skill_nmme = 1 - [0.15890   0.18053   0.18078   0.18216   0.18248   0.18308   0.18512]/0.26274;
p = plot([0 6], [skill_trend skill_trend], ':', 0:6, skill_nmme);
set (gca, 'linewidth', lw);
axis([0 6 0 0.4])
title("Global land")
ylabel('Normalized log likelihood')
set (gca,'position',[0.15 0.2 0.8 0.7]);
l = legend ('Trend', 'NMME', "location", "southwest");
set (gca, 'fontsize', fs);
text (texx, texy, ["(" letter ")"], "fontsize", fs, "units", "normalized")
print(['-d' plot_format], [plot_dir 'gl_ig' '.' plot_suffix])

#NA
clf('reset');
letter = 'b';
skill_trend = 1 - 0.15764/0.18136;
skill_nmme = 1 - [0.12458   0.14461   0.14645   0.14933   0.14844   0.15356   0.15117]/0.18136;
p = plot([0 6], [skill_trend skill_trend], ':', 0:6, skill_nmme);
set (gca, 'linewidth', lw);
axis([0 6 0 0.4])
title("North America")
xlabel("Forecast lead (months)")
ylabel('Normalized log likelihood')
set (gca,'position',[0.15 0.2 0.8 0.7]);
#l = legend ('Trend', 'NMME', "location", "southwest");
set (gca, 'fontsize', fs);
text (texx, texy, ["(" letter ")"], "fontsize", fs, "units", "normalized")
print(['-d' plot_format], [plot_dir 'na_ig' '.' plot_suffix])


