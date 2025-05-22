## Copyright (C) 2016 Nir Krakauer
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

##Run analyses and make plots for "Temperature trends and prediction skill in NMME seasonal forecasts"
##
##Note: data_dir, out_dir, plot_dir should be defined as variables in the current workspace before running this script

## Author: Nir Krakauer <nkrakauer@ccny.cuny.edu>

pkg load nan splines netcdf statistics

  prefix_dir = [home_dir];
  data_dir = [prefix_dir filesep() 'other/data/NMME/'];   
  out_dir = [data_dir "out" filesep()];
  plot_dir = [data_dir "plots" filesep()];

if !exist("data_dir", "var")
  warning ("data_dir not defined, using current directory")
  data_dir = [pwd filesep()];
endif
if !exist("out_dir", "var")
  out_dir = [data_dir "out" filesep()];
  warning (["out_dir not defined, using " out_dir])
endif

re_run = false;
draw_figures = false;


nl = 12; #lags 0.5 to 11.5 [0 to 11]
gcms = {'CMC1-CanCM3', 'CMC2-CanCM4', 'COLA-RSMAS-CCSM3', 'COLA-RSMAS-CCSM4', 'GFDL-CM2p1-aer04', ...
        'GFDL-CM2p5-FLOR-A06', 'GFDL-CM2p5-FLOR-B01', 'NASA-GMAO-062012', 'NCEP-CFSv2'};
n_gcms = numel(gcms);
obs_dataset = 'best';

subsets = {'land', 'ocean'}
n_subsets = numel (subsets);

years = 2012:2015;
months = 1:12;

opts = struct ('data_dir', data_dir, ...
                  're_download', false, ...
                  're_process', false, ...
                  'var_name', 'tref', ...
                  'fcst_lons', (0:359)', ...
                  'fcst_lats', (-90:90)', ... 
                  'fcst_start_year', 1982, ...
                  'obs_start_year', 1957, ...
                  'obs_end_year', 2015, ...                  
                  'obs_dataset', obs_dataset, ...
                  'p_spline', 1E-4, ...
                  'dof_spline', 2.5, ...
                  'lamda', 10, ...
                  'predict_adj', false, ...
                  'predict_adj_years', 1995:1999, ...
                  'predict_adj_months', 1:12, ...
                  'clim_years', 1981:2010, ...
                  'gcms', {gcms} ...
                  );


for l = 1:nl
  lag = l - 1;

  for s = 1:n_subsets

    subset = subsets{s};

    file_out = [out_dir 'nmme_trend_results'  '_' num2str(lag, '%2.2i') '_' obs_dataset '_' subset '.mat']; #e.g. nmme_trend_results_00_best_land.mat

    opts.subset_verify = subset;

    if !exist(file_out, "file") || re_run

      #compute trends in models and observations, and mean correlations between them
      start_year = opts.fcst_start_year; #1982
      n_years = 34;
      end_year = start_year + n_years - 1;
      n_perms = 2;

      #observations
      m = 360*181; #grid cells
      obs_all = nan (m, 12*n_years, "single");
      for month = 1:12
        file_obs = sefo_obs_assemble (end_year+1, month, opts);
        load (file_obs, "obs")
        obs_all(:, month:12:end) = obs(:, (end-n_years+1):end);
      endfor
      disp(["assembled obs, lag = " num2str(lag) ", s = " num2str(s)])
      obs = obs_all; clear obs_all
      re_download = opts.re_download;
      fcst_lats = opts.fcst_lats;
      fcst_lons = opts.fcst_lons;
      nlats = numel(fcst_lats);
      nlons = numel(fcst_lons);
      fcst_lon_lims = nan(nlons+1, 1);
      fcst_lon_lims(2:(end-1)) = (fcst_lons(1:(end-1)) + fcst_lons(2:end))/2;
      fcst_lon_lims(1) = (fcst_lons(1) + fcst_lons(end))/2 - 180;
      fcst_lon_lims(end) = fcst_lon_lims(1) + 360;
      fcst_lat_lims = nan(nlats+1, 1);
      fcst_lat_lims(2:(end-1)) = (fcst_lats(1:(end-1)) + fcst_lats(2:end))/2;
      fcst_lat_lims(1) = -90;
      fcst_lat_lims(end) = 90;
      fcst_lat_lims = sind (fcst_lat_lims);
      #get a land mask from BEST and regrid it to match forecast (NMME)
      source_file_name = ['Land_and_Ocean_LatLong1.nc'];
      if ~exist([data_dir source_file_name], 'file') || re_download
        data_url = ['http://berkeleyearth.lbl.gov/auto/Global/Gridded/' source_file_name];
        urlwrite (data_url, [data_dir source_file_name])  #download and save netcdf file
      endif
      source_file = [data_dir source_file_name];
      lon_obs = ncread (source_file, 'longitude'); #-179.5:179.5
      lat_obs = ncread (source_file, 'latitude'); #-89.5:89.5
      land_mask = ncread (source_file, 'land_mask');
      arr = land_mask';
      arr = [arr arr arr];
      obs_lat_lims = sind([-90 lat_obs(1:(end-1))'+diff(lat_obs')/2 90])';
      obs_lon_lims = [lon_obs(end)-360 lon_obs' lon_obs(1)+360];
      obs_lon_lims =  obs_lon_lims(1:(end-1))'+diff(obs_lon_lims')/2;
      obs_lon_lims = [obs_lon_lims(1:(end-1))'-360 obs_lon_lims' obs_lon_lims(2:end)'+360]';
      coords = {obs_lat_lims, obs_lon_lims};
      new_coords = {fcst_lat_lims, fcst_lon_lims};
      land_mask_regrid = regrid (coords, arr, new_coords);
      area_wgts = diff(fcst_lat_lims) * ones(1, nlons);
      #subset land or ocean
      switch subset
        case 'land'
          w = area_wgts(:) .* land_mask_regrid(:);
        case 'ocean'
          w = area_wgts(:) .* (1 - land_mask_regrid(:));
        otherwise
          error (['Undefined subset: ' subset])
      endswitch
      obs_g = mean (double(obs), 1, double(w));
      X_m = zeros(12*n_years, 12);
      for i = 1:12
        X_m(i:12:end, i) = 1;
      endfor
      X_t = (1:(12*n_years))';
      X_t -= mean(X_t);
      y = obs_g';
      #[b_obs, bint_obs, r_obs, ~, stats_obs] = regress(y(2:end), [X_m(2:end, :) X_t(2:end) y(1:(end-1))]); #better model, but slope is harder to interpret
      disp(["begin regressions, lag = " num2str(lag) ", s = " num2str(s)])
      [b_obs, bint_obs, r_obs, ~, stats_obs] = regress(y, [X_m X_t]);
      B = reshape (double(obs), [m*12 n_years])';
      B_detrend = detrend (B, 1);
      #forecasts
      b_fcst = nan (13, n_gcms+1);
      bint_fcst = nan (13, 2, n_gcms+1);
      corr_mean = corr_mean_detrend = nan (n_gcms+1, 1);
      corr_mean_p = corr_mean_detrend_p = nan (n_gcms+1, n_perms);
      fcsts = nan (m, 12*n_years, n_gcms, "single");
      for month = 1:12
        file = sefo_fcst_assemble (end_year, month, lag, opts);
        load (file, "fcsts_mean", "nrs", "missing_models_excluded")
        fcsts(:, month:12:end, ~missing_models_excluded) = fcsts_mean;
      endfor
      clear fcsts_mean
      fcsts(:, :, n_gcms+1) = mean (fcsts(:, :, 1:n_gcms), 3); #the multimodel mean
      fcsts_g_uw = mean (fcsts, 1); #unweighted mean of grid cells
      for month = 1:12
        month
        disp (squeeze(fcsts_g_uw)(month:12:end, :))
      endfor
      fcsts_g = mean (fcsts, 1, double(w));
      fcsts_g_deseas = squeeze(fcsts_g) - reshape(repmat(mean(reshape(squeeze(fcsts_g), 12, n_years, n_gcms+1), 2), [1 n_years 1]), 12*n_years, n_gcms+1);
      for i = 1:(n_gcms+1)
        y = squeeze(fcsts_g)(:, i);
        if all(isfinite(y))
          [b, bint, ~, ~, stats] = regress(y, [X_m X_t]);
          b_fcst(:, i) = b;
          bint_fcst(:, :, i) = bint;
        elseif i <= n_gcms
          warning (['missing data for model ' gcms{i} ', lag ' num2str(lag+0.5)])
        endif

        #correlation between models and observations
        A = reshape(double(fcsts(:, :, i)), [m*12 n_years])';
        corr_mean(i) = mean(corrc (A, B), 2, double(w*ones(1, 12, "single"))(:));

        #after removal of linear trend
        A_detrend = detrend (A, 1);
        corr_mean_detrend(i) = mean(corrc (A_detrend, B_detrend), 2, double(w*ones(1, 12, "single"))(:));

        for p = 1:n_perms
          A_p = A(randperm(n_years), :);
          corr_mean_p(i, p) = mean(corrc (A_p, B), 2, double(w*ones(1, 12, "single"))(:));

          A_p_detrend = detrend (A_p, 1);
          corr_mean_p_detrend(i, p) = mean(corrc (A_p_detrend, B_detrend), 2, double(w*ones(1, 12, "single"))(:));
        endfor

      endfor
      clear fcsts
      
      save (file_out, 'n_years', 'n_gcms', 'obs_g', 'b_obs', 'bint_obs', 'r_obs', 'stats_obs', 'b_fcst', 'bint_fcst', 'gcms', 'corr_mean', 'corr_mean_detrend', 'corr_mean_p', 'corr_mean_p_detrend', 'fcst_lons', 'fcst_lats', 'fcsts_g_uw', 'fcsts_g', 'fcsts_g_deseas')
    
    endif


    start_year = opts.fcst_start_year; #1982
    n_years = 34;
    end_year = start_year + n_years - 1;


    if lag == 0
      tag = ['_' obs_dataset '_' subset];
      file_out = [out_dir 'sefo_verify_script_out_hist' tag '.mat'];


      if !exist(file_out, "file") || re_run
        tic;

        #verify forecast methods that are only based on historical observations
        system(['rm ' data_dir 'verify*.mat']);
        methods = {'C', 'MA', 'EW', 'EW-a'};
        #MA, T, EW, EW-a
        opts.methods = methods;
        file = sefo_verify (years, months, lag, opts);
        load (file)

        t_taken = toc;
        disp (t_taken)

        save (file_out, "years", "months", "lag", "methods", "rmse", "bias", "nll", "crps", "predent", "ks", "hists", "factors_adj", "rmse_map", "bias_map", "nll_map", "predent_map", "rmse_ts", "bias_ts", "nll_ts", "predent_ts", "nll_cats", "prob_cats", "freq_obs", "freq_for", "ll_for")
      endif
    endif

    tic;


    tag = ['_' obs_dataset '_' subset '_' num2str(lag, '%2.2i')];
    file_out = [out_dir 'sefo_verify_script_out_nmme' tag '.mat'];
    if !exist(file_out, "file") || re_run
      #methods that employ NMME outputs
      system(['rm ' data_dir 'verify*.mat']);
      methods = {'M', 'MT', 'MS', 'MST', 'MM', 'MMT'}; #M = model, T = (linear) trend, MM = multimodel, S = scaled
                            #M, MT, MS, MST, MM, MMT
      opts.methods = methods;
      file = sefo_verify (years, months, lag, opts);
      load (file)

      t_taken = toc;
      disp (t_taken)

      save (file_out, "years", "months", "lag", "methods", "rmse", "bias", "nll", "crps", "predent", "ks", "hists", "factors_adj", "rmse_map", "bias_map", "nll_map", "predent_map", "rmse_ts", "bias_ts", "nll_ts", "predent_ts", "nll_cats", "prob_cats", "freq_obs", "freq_for", "ll_for")
    endif
  endfor
  
  #remove all intermediate files
  strng = ['for i in ' data_dir '*.*; do rm "$i"; done'];
  ##system (strng);
  
endfor

if draw_figures

  if !exist("plot_dir", "var")
    warning ("plot_dir not defined, using current directory")
    plot_dir = [pwd filesep()];;
  endif

  obs_dataset = 'best';
  subsets = {'land', 'ocean'}
  n_subsets = numel (subsets);
  nl = 12; #0.5 to 11.5 [0 to 11]
  b_fcsts = bint_fcsts = corr_means = corr_mean_detrends = b_obss = bint_obss = fcstss_g = fcstss_g_uw = fcstss_g_deseas = obss_g_deseas = ...
    rmses = biases = predents = nlls = crpss = kss = nll_maps = rmse_maps = bias_maps = predent_maps = nll_tss = rmse_tss = histss = nll_catss = prob_catss = freq_fors = ll_fors = freq_obss = ...
    rmsesm = biasesm = predentsm = nllsm = crpsm = kssm = nll_mapsm = rmse_mapsm = bias_mapsm = predent_mapsm = nll_tssm = rmse_tssm = histssm = nll_catssm = prob_catssm = freq_forsm = ll_forsm = corrs_mean_p = corrs_mean_p_detrend = ...
    [];
  for i = 1:nl
    lag = i - 1;
    for s = 1:n_subsets
      subset = subsets{s};
      file_out = [out_dir 'nmme_trend_results'  '_' num2str(lag, '%2.2i') '_' obs_dataset '_' subset '.mat'];
      load (file_out)
      b_fcsts(:, :, i, s) = b_fcst;
      bint_fcsts(:, :, :, i, s) = bint_fcst;
      corr_means(:, i, s) = corr_mean;
      corr_mean_detrends(:, i, s) = corr_mean_detrend;
      fcstss_g(:, :, i, s) = fcsts_g;
      fcstss_g_uw(:, :, i, s) = fcsts_g_uw;
      fcstss_g_deseas(:, :, i, s) = fcsts_g_deseas;
      corrs_mean_p(:, :, i, s) = corr_mean_p;
      corrs_mean_p_detrend(:, :, i, s) = corr_mean_p_detrend;
      if i == 1
        b_obss(:, s) = b_obs;
        bint_obss(:, :, s) = bint_obs;
        obs_g_deseas = obs_g' - reshape(repmat(mean(reshape(obs_g, 12, n_years), 2), [1 n_years]), 12*n_years, 1);
        obss_g_deseas(:, s) = obs_g_deseas;

        tag = ['_' obs_dataset '_' subset];
        file_out = [out_dir 'sefo_verify_script_out_hist' tag '.mat'];
        load (file_out)
        rmses(:, s) = rmse; biases(:, s) = bias; predents(:, s) = predent;
        nlls(:, s) = nll; crpss(:, s) = crps; kss(:, s) = ks;
        nll_maps(:, :, s) = nll_map;
        rmse_maps(:, :, s) = rmse_map;
        bias_maps(:, :, s) = bias_map;
        predent_maps(:, :, s) = predent_map;
        nll_tss(:, :, s) = nll_ts;
        rmse_tss(:, :, s) = rmse_ts;
        histss(:, :, s) = hists;
        nll_catss(:, :, s) = nll_cats;
        prob_catss(:, :, s) = prob_cats;
        freq_fors(:, :, s) = freq_for;
        ll_fors(:, :, s) = ll_for;
        freq_obss(:, s) = freq_obs;
      endif
    
      tag = ['_' obs_dataset '_' subset '_' num2str(lag, '%2.2i')];
      file_out = [out_dir 'sefo_verify_script_out_nmme' tag '.mat'];
      load (file_out)
      rmsesm(:, i, s) = rmse; biasesm(:, i, s) = bias; predentsm(:, i, s) = predent;
      nllsm(:, i, s) = nll; crpssm(:, i, s) = crps; kssm(:, i, s) = ks;
      nll_mapsm(:, :, i, s) = nll_map;
      rmse_mapsm(:, :, i, s) = rmse_map;
      bias_mapsm(:, :, i, s) = bias_map;
      predent_mapsm(:, :, i, s) = predent_map;
      nll_tssm(:, :, i, s) = nll_ts;
      rmse_tssm(:, :, i, s) = rmse_ts;
      histssm(:, :, i, s) = hists;
      nll_catssm(:, :, i, s) = nll_cats;
      prob_catssm(:, :, i, s) = prob_cats;
      freq_forsm(:, :, i, s) = freq_for;
      ll_forsm(:, :, i, s) = ll_for;
    endfor
  endfor
    
    lw = 2; %line width
    fnt = 15; %font size
    fnt_l = 12; %font size for legend
    set (0,"Defaulttextfontsize", fnt); 
    set(0,"Defaultaxesfontsize", fnt);
    set (0, "defaultlinelinewidth", lw);

    y = 120*[squeeze(b_fcsts(13, :, :, 1)); b_obss(13, 1)*ones(1, 12)];
    p = plot (0.5:11.5, y);
    axis ([0.25 11.75 min(y(:))-0.01 max(y(:))+0.01])
    set (p(end-1), 'linewidth', 2*lw)
    set (p(end-1), 'color', [0 0 1])
    set (p(end), 'linewidth', 2*lw)
    set (p(end), 'color', [0 0 0])
    set (p(1:3:(end-1)), "linestyle", "--")
    set (p(2:3:(end-1)), "linestyle", ":")
    set (p(3:3:(end-1)), "linestyle", "-.")
    l = legend ('CMC1-CanCM3', 'CMC2-CanCM4', 'COLA-RSMAS-CCSM3', 'COLA-RSMAS-CCSM4', 'GFDL-CM2p1-aer04', 'GFDL-CM2p5-FLOR-A06', 'GFDL-CM2p5-FLOR-B01', 'NASA-GMAO-062012', 'NCEP-CFSv2', 'NMME-avg', 'BEST', "location", "eastoutside");
    set (l, 'fontsize', fnt_l);
    legend boxoff
    xlabel ('Prediction lag (months)')
    ylabel ('Warming trend (K/decade)')
    print ('-depsc', '-Fcourier', [plot_dir 'trend'])


    lag_index = 2;
    subset_index = 1;
    A = [fcstss_g_deseas(:, 1:n_gcms, lag_index, subset_index) obss_g_deseas(:, subset_index)+1];
    p = plot (1982+((1:(12*n_years))-0.5)/12, A);
    axis ([1981.5 2016.5 min(A(:))-0.1 max(A(:))+0.1])
    set (p(end), 'color', [0 0 0])
    set (p(1:3:(end-1)), "linestyle", "--")
    set (p(2:3:(end-1)), "linestyle", ":")
    set (p(3:3:(end-1)), "linestyle", "-.")
    ylabel ('Temperature anomaly (K)')
    grid on
    l = legend ('CMC1-CanCM3', 'CMC2-CanCM4', 'COLA-RSMAS-CCSM3', 'COLA-RSMAS-CCSM4', 'GFDL-CM2p1-aer04', 'GFDL-CM2p5-FLOR-A06', 'GFDL-CM2p5-FLOR-B01', 'NASA-GMAO-062012', 'NCEP-CFSv2', 'obs', "location", "eastoutside");
    set (l, 'fontsize', fnt_l);
    legend boxoff
    print ('-depsc', '-Fcourier', [plot_dir 'T_ts'])

    p = plot (0.5:11.5, corr_means(:, :, 1));
    axis ([0.25 11.75 0 0.6])
    set (p(end), 'linewidth', 2*lw)
    set (p(end), 'color', [0 0 1])
    set (p(1:3:(end-1)), "linestyle", "--")
    set (p(2:3:(end-1)), "linestyle", ":")
    set (p(3:3:(end-1)), "linestyle", "-.")
    l = legend ('CMC1-CanCM3', 'CMC2-CanCM4', 'COLA-RSMAS-CCSM3', 'COLA-RSMAS-CCSM4', 'GFDL-CM2p1-aer04', 'GFDL-CM2p5-FLOR-A06', 'GFDL-CM2p5-FLOR-B01', 'NASA-GMAO-062012', 'NCEP-CFSv2', 'NMME-avg', "location", "northeast");
    set (l, 'fontsize', fnt_l);
    legend boxoff
    xlabel ('Prediction lag (months)')
    ylabel ('Mean correlation')
    print('-depsc', '-Fcourier', [plot_dir 'corr'])

    p = plot (0.5:11.5, corr_mean_detrends(:, :, 1));
    axis ([0.25 11.75 0 0.6])
    set (p(end), 'linewidth', 2*lw)
    set (p(end), 'color', [0 0 1])
    set (p(1:3:(end-1)), "linestyle", "--")
    set (p(2:3:(end-1)), "linestyle", ":")
    set (p(3:3:(end-1)), "linestyle", "-.")
    xlabel ('Prediction lag (months)')
    ylabel ('Mean correlation')
    print('-depsc', '-Fcourier', [plot_dir 'corr_detrend'])
    
    for s = 1:2
      if s == 1
        tag = '';
      else
        tag = '_ocean';
      endif
      A = [nllsm(1:6, :, s); nlls(4, s)*ones(1, nl)];
      p = plot (0.5:11.5, A);
      axis ([0.25 11.75 min(A(:))-0.01 max(A(:))+0.01])
      set (p(end), 'linewidth', 2*lw)
      set (p(end), 'color', [0 0 0])
      set (p(1:3:(end-1)), "linestyle", "--")
      set (p(2:3:(end-1)), "linestyle", ":")
      set (p(3:3:(end-1)), "linestyle", "-.")
      l = legend ('M', 'MT', 'MS', 'MST', 'MM', 'MMT', 'EW-a', "location", "southeast");
      set (l, 'fontsize', fnt_l);
      legend boxoff
      xlabel ('Prediction lag (months)')
      ylabel ('NLL (nats)')
      print('-depsc', '-Fcourier', [plot_dir 'nll' tag])

      A = [rmsesm(1:6, :, s); rmses(4, s)*ones(1, nl)];
      p = plot (0.5:11.5, A);
      axis ([0.25 11.75 min(A(:))-0.01 max(A(:))+0.01])
      set (p(end), 'linewidth', 2*lw)
      set (p(end), 'color', [0 0 0])
      set (p(1:3:(end-1)), "linestyle", "--")
      set (p(2:3:(end-1)), "linestyle", ":")
      set (p(3:3:(end-1)), "linestyle", "-.")
      l = legend ('M', 'MT', 'MS', 'MST', 'MM', 'MMT', 'EW-a', "location", "southeast");
      set (l, 'fontsize', fnt_l);
      legend boxoff
      xlabel ('Prediction lag (months)')
      ylabel ('RMSE (K)')
      print('-depsc', '-Fcourier', [plot_dir 'rmse' tag])
    endfor
    
    A = 100*[freq_fors(:, 1:4, 1) freq_obss(:, 1)];
    p = plot (0.5:99.5, A);
    axis ([0 100 0 max(A(:))+0.2])
    set (p(end), 'linewidth', 2*lw)
    set (p(end), 'color', [0 0 0])
    set (p(1:3:(end-1)), "linestyle", "--")
    set (p(2:3:(end-1)), "linestyle", ":")
    set (p(3:3:(end-1)), "linestyle", "-.")
    l = legend ('C', 'MA', 'EW', 'EW-a', 'obs', "location", "northwest");
    set (l, 'fontsize', fnt_l);
    legend boxoff
    xlabel ('Climatology percentile')
    ylabel ('Mean probability (%)')
    print('-depsc', '-Fcourier', [plot_dir 'freq_obs'])

    A = 100*[freq_forsm(:, 1:6, 2, 1) freq_obss(:, 1)];
    p = plot (0.5:99.5, A);
    axis ([0 100 0 max(A(:))+0.2])
    set (p(end), 'linewidth', 2*lw)
    set (p(end), 'color', [0 0 0])
    set (p(1:3:(end-1)), "linestyle", "--")
    set (p(2:3:(end-1)), "linestyle", ":")
    set (p(3:3:(end-1)), "linestyle", "-.")
    l = legend ('M', 'MT', 'MS', 'MST', 'MM', 'MMT', 'obs', "location", "northwest");
    set (l, 'fontsize', fnt_l);
    legend boxoff
    xlabel ('Climatology percentile')
    ylabel ('Mean probability (%)')
    print('-depsc', '-Fcourier', [plot_dir 'freq_nmme'])

    A = 100*[freq_forsm(:, 6, 2, 1) prob_catssm(:, 6, 2, 1)];
    p = plot (0.5:99.5, A);
    axis ([0 100 0 max(A(:))+0.2])
    l = legend ('unconditional', 'conditional', "location", "northwest");
    set (l, 'fontsize', fnt_l);
    legend boxoff
    xlabel ('Climatology percentile')
    ylabel ('Mean probability (\%)')
    print('-depsc', '-Fcourier', [plot_dir 'cond_mmt'])

    A = [ll_forsm(:, 6, 2, 1) -nll_catssm(:, 6, 2, 1)];
    #A = [ll_fors(:, 4, 1) -nll_catss(:, 4, 1)];
    p = plot (0.5:99.5, A);
    axis ([0 100 min(A(:))-0.1 max(A(:))+0.1])
    l = legend ('unconditional', 'conditional', "location", "northwest");
    #set (l, 'fontsize', fnt_l);
    legend boxoff
    xlabel ('Climatology percentile')
    ylabel ('Mean log(p)')
    print('-depsc', '-Fcourier', [plot_dir 'cond_ll_mmt'])

    #IG from MMT vs. EW-a
    A = nll_maps(:, 4, 1) - nll_mapsm(:, 6, 2, 1);
    colormap("default"); imagesc (fcst_lons, fcst_lats, reshape(A, 181, 360), [-0.5 0.5]); axis xy; pcontinents2; colorbar
    print('-depsc', '-Fcourier', [plot_dir 'map_ig'])

endif

%{
#find mean posterior versus prior NLL for hot/cold relative extremes for different absolute climatological temperatures
#MMT, lag 1.5

pkg load netcdf
nl = 12; #lags 0.5 to 11.5 [0 to 11]
gcms = {'CMC1-CanCM3', 'CMC2-CanCM4', 'COLA-RSMAS-CCSM3', 'COLA-RSMAS-CCSM4', 'GFDL-CM2p1-aer04', 'GFDL-CM2p5-FLOR-A06', 'GFDL-CM2p5-FLOR-B01', 'NASA-GMAO-062012', 'NCEP-CFSv2'};
n_gcms = numel(gcms);
obs_dataset = 'best';

years = 2012:2015;
months = 1:12;
lag = 1;

opts = struct ('data_dir', data_dir, ...
                  're_download', false, ...
                  're_process', false, ...
                  'var_name', 'tref', ...
                  'fcst_lons', (0:359)', ...
                  'fcst_lats', (-90:90)', ... 
                  'fcst_start_year', 1982, ...
                  'obs_start_year', 1957, ...
                  'obs_end_year', 2015, ...                  
                  'obs_dataset', obs_dataset, ...
                  'p_spline', 1E-4, ...
                  'dof_spline', 2.5, ...
                  'lamda', 10, ...
                  'predict_adj', false, ...
                  'clim_years', 1981:2010, ...
                  'gcms', {gcms}, ...
                  'methods', {{'MMT'}}, ...
                  'subset_verify', '' ...
                  );

source_file_name = ['Land_and_Ocean_LatLong1.nc'];
source_file = [data_dir source_file_name];
lon_obs = ncread (source_file, 'longitude');
lat_obs = ncread (source_file, 'latitude');
clim_best = ncread (source_file, 'climatology'); #360x180x12

clim_nbins = 100;

fcst_lons = getfield (opts, 'fcst_lons');
fcst_lats = getfield (opts, 'fcst_lats');
nlats = numel(fcst_lats);
nlons = numel(fcst_lons);
fcst_lon_lims = nan(nlons+1, 1);
fcst_lon_lims(2:(end-1)) = (fcst_lons(1:(end-1)) + fcst_lons(2:end))/2;
fcst_lon_lims(1) = (fcst_lons(1) + fcst_lons(end))/2 - 180;
fcst_lon_lims(end) = fcst_lon_lims(1) + 360;
fcst_lat_lims = nan(nlats+1, 1);
fcst_lat_lims(2:(end-1)) = (fcst_lats(1:(end-1)) + fcst_lats(2:end))/2;
fcst_lat_lims(1) = -90;
fcst_lat_lims(end) = 90;
fcst_lat_lims = sind (fcst_lat_lims);
obs_lat_lims = sind([-90 lat_obs(1:(end-1))'+diff(lat_obs')/2 90])';
obs_lon_lims = [lon_obs(end)-360 lon_obs' lon_obs(1)+360];
obs_lon_lims =  obs_lon_lims(1:(end-1))'+diff(obs_lon_lims')/2;
obs_lon_lims = [obs_lon_lims(1:(end-1))'-360 obs_lon_lims' obs_lon_lims(2:end)'+360]';
coords = {obs_lat_lims, obs_lon_lims};
new_coords = {fcst_lat_lims, fcst_lon_lims};
arr_regrid = [];
for m = 1:12
  arr = clim_best(:, :, m)';
  arr_aug = [arr arr arr];
  arr_regrid(:, :, m) = regrid (coords, arr_aug, new_coords);
endfor
clim_best = arr_regrid;

land_mask = ncread (source_file, 'land_mask');
arr_regrid = [];
for m = 1:12
  arr = land_mask';
  arr_aug = [arr arr arr];
  arr_regrid(:, :, m) = regrid (coords, arr_aug, new_coords);
endfor
land_mask = arr_regrid;


clim_q = quantile (clim_best(land_mask > 0.5), (1:(clim_nbins-1))/clim_nbins);
clim_p = reshape (interp1(clim_q, 2:clim_nbins, clim_best(:), "previous", "extrap"), size(clim_best));
clim_p(clim_best < clim_q(1)) = 1;
clim_q_mean = accumarray (clim_p(:), clim_best(:), [], @mean);



#obs_cats, thresh_probs, prob_cat

#commands based on sefo_verify
#compile observations
clim_years = getfield (opts, 'clim_years');
rank_calc = true;
n_years = numel(years);
n_months = numel(months);
years = ones(n_months, 1) * years(:)';
months = months(:) * ones(1, n_years);
nt = n_years*n_months;
n_months_combine = 1;

obss = area_wgtss = [];
for t = 1:nt
  year = years(t);
  month = months(t);
  predict_month = month + (0:(n_months_combine-1))';
  predict_year = year * ones(n_months_combine, 1);
  inds_next_year = find(predict_month > 12);
  predict_year(inds_next_year) = predict_year(inds_next_year) + 1;
  predict_month(inds_next_year) = predict_month(inds_next_year) - 12;
  n_days_per_month = days_per_month (predict_month, predict_year);
  w = n_days_per_month / sum(n_days_per_month);
  obs = 0;
  for m = 1:n_months_combine
    file_obs = sefo_obs_read (predict_year(m), predict_month(m), opts);
    load (file_obs)
    if ~exist("obs_inds", "var") || isempty(obs_inds)
      obs_inds = isfinite (arr_regrid);
    endif
    obs += w(m)*arr_regrid(obs_inds);
    clear arr_regrid
  endfor
  if t == 1
    np = numel (obs);
    obss = nan (np, nt, "single");
  endif
  obss(:, t) = obs;
  
  if rank_calc
  
    #get climatology
    clim_opts = opts;
    clim_opts.obs_start_year = min(clim_years);
    clim_obs_file = sefo_obs_assemble (max(clim_years)+1, month, clim_opts);
    load (clim_obs_file, "obs");
    
    #based on assumed normality over the climatology period
    n_cats = 100;
    if t == 1
      obs_cats = zeros (np, nt, "uint16");
      thresh_obs_cats = nan (np, nt, 2, "single");
      thresh_cats = nan (np, nt, n_cats-1, "single");
    endif
    thresh_z = tinv ((1:(n_cats-1))'/n_cats, n_cats-1);
    clim_mean = mean(obs, 2);
    clim_sd = std(obs, [], 2);
    obs_cat = ceil (n_cats * tcdf ((obss(:, t) - clim_mean) ./ clim_sd, n_cats-1));
    obs_cat(obs_cat == 0) = 1; #cases where the z score is low enough that tcdf underflows to 0
    thresh_obs_cat = nan (np, 2);
    inds = (obs_cat == 1);
    thresh_obs_cat(!inds, 1) = clim_mean(!inds) + thresh_z(obs_cat(!inds) - 1) .* clim_sd(!inds);
    thresh_obs_cat(inds, 1) = -Inf;
    inds = (obs_cat == n_cats);
    thresh_obs_cat(!inds, 2) = clim_mean(!inds) + thresh_z(obs_cat(!inds)) .* clim_sd(!inds);
    thresh_obs_cat(inds, 2) = Inf;
            
    obs_cats(:, t) = obs_cat;
    thresh_obs_cats(:, t, :) = thresh_obs_cat;
    thresh_cats(:, t, :) = clim_mean + clim_sd * thresh_z';
    disp ([datestr(clock) "\n" 'computed observation ranks'])
  endif
  disp ([datestr(clock) "\n" 'compiled observations, t = ' num2str(t)])
endfor

method = getfield(opts, 'methods'){1};
for t = 1:nt
  year = years(t);
  month = months(t);
  predict_month = month + (0:(n_months_combine-1))';
  predict_year = year * ones(n_months_combine, 1);
  inds_next_year = find(predict_month > 12);
  predict_year(inds_next_year) = predict_year(inds_next_year) + 1;
  predict_month(inds_next_year) = predict_month(inds_next_year) - 12;
  file_predict = sefo_predict (predict_year, predict_month, lag, method, opts){1};
  load (file_predict, "mu", "sigma", "nu", "area_wgts")
  predict_mu(:, t) = mu(:);
  predict_sigma(:, t) = sigma(:);
  predict_nu(:, t) = nu(:);
  if m == 1
    area_wgts = area_wgts / mean(double(area_wgts));
    area_wgtss(:, t) = area_wgts;
  endif
  disp ([datestr(clock) "\n" 'compiled predictions, t = ' num2str(t)])
endfor
disp ([datestr(clock) "\n" 'collected predictions'])

#forecast probability of each climatology category
w = mean (double (area_wgtss), 2);
#if rank_calc
  #thresh_probs = diff (cat (3, zeros(np, nt, 1), tcdf ((thresh_cats - predict_mu) ./ predict_sigma, repmat(single(predict_nu), [np 1 n_cats-1])), ones(np, nt, 1)), [], 3);
  #freq_for = mean (mean (thresh_probs, 1, double(w)), 2);
  #ll_for = mean (mean (log(thresh_probs), 1, double(w)), 2);
#endif


clim_p = reshape (clim_p, np, 12);
land_mask = reshape (land_mask, np, 12);

l_min = exp (-40); #a small likelihood threshold for eliminating zero-probability values that are due to roundoff
ll_uncond_clim = ll_cond_clim = n_cond_clim = nan (clim_nbins, n_cats);
n_cond_clim = zeros (clim_nbins, n_cats);
min_thresh = nan (clim_nbins, 1);
for i = 1:clim_nbins
  inds = repmat((clim_p == i) & (land_mask > 0.5), 1, n_years);
  A = reshape(thresh_cats, np*nt, n_cats-1)(inds, :);
  ni = size (A, 1);
  B = (A - predict_mu(inds)) ./ predict_sigma(inds);
  C = reshape(repmat(single(predict_nu), [np 1 n_cats-1]), np*nt, n_cats-1)(inds, :);
  thresh_probs = diff (cat (2, zeros(ni, 1), tcdf (double(B), C), ones(ni, 1)), [], 2);
  min_thresh(i) = min (thresh_probs(:));
  thresh_probs = max (thresh_probs, l_min);
  ll_uncond_clim(i, :) = mean(log(thresh_probs), 1);
  cat_s = obs_cats(inds);
  for j = 1:n_cats
    if any (cat_s == j)
      ll_cond_clim(i, j) = mean (log(thresh_probs(cat_s == j, j)));
      n_cond_clim(i, j) = sum (cat_s == j);
    endif
  endfor
endfor

    lw = 2; %line width
    fnt = 15; %font size
    fnt_l = 12; %font size for legend
    set (0,"Defaulttextfontsize", fnt); 
    set(0,"Defaultaxesfontsize", fnt);
    set (0, "Defaultlinelinewidth", lw);

    a(:, 1) = mean(ll_cond_clim(1:50, :) - ll_uncond_clim(1:50, :))';
    a(:, 2) = mean(ll_cond_clim(51:87, :) - ll_uncond_clim(51:87, :))';
    a(:, 3) = mean(ll_cond_clim(88:end, :) - ll_uncond_clim(88:end, :))';
    p = plot (0.5:99.5, a);
    axis ([0 100 0 max(a(:))+0.1])
    l = legend ('cold', 'temperate', 'hot', "location", "north");
    set (l, 'fontsize', fnt_l);
    set (p, 'linewidth', lw);
    set (p(2), 'color', [0 1 0]) #green
    legend boxoff
    xlabel ('Climatology percentile')
    ylabel ('Mean log(p) change')
    print('-depsc', '-Fcourier', [plot_dir 'll_t_mmt'])

%}





%{
    #observation rank histogram [not area-weighted]
    plot([histss(:, 4, 1) histssm(:, 6, 2, 1)]*(100/(65160*48)))
  
    #mean by latutude
    plot(fcst_lats, mean(reshape(A, 181, 360)'))
    #models best near Equator, worst in Antarctica

%}
